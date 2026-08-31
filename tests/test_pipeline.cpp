// Exercises every stage except the BAM reader with synthetic records, so the
// classification / detection / EM / output path can be tested without htslib.
//   g++ -std=c++17 -O2 -I../src test_pipeline.cpp -o test_pipeline && ./test_pipeline
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "common.hpp"
#include "detect.hpp"
#include "group.hpp"
#include "gtf.hpp"
#include "output.hpp"
#include "pipeline.hpp"
#include "quant.hpp"
#include "record.hpp"

using namespace brocoli;

static int g_fail = 0;
#define CHECK(cond, msg)                                                       \
    do {                                                                       \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++g_fail; }       \
        else std::cout << "  ok: " << msg << "\n";                             \
    } while (0)

static const char* kGtf = R"(chr1	src	exon	100	200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	src	exon	300	400	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	src	exon	500	600	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	src	exon	100	200	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	src	exon	500	600	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	src	exon	1000	1100	.	+	.	gene_id "G2"; transcript_id "T3";
)";

static void addRead(std::vector<ReadRec>& v, const std::string& name,
                    int beg, int end, const IntervalVec& sjs, int file = 0,
                    const std::string& bc = "", const std::string& umi = "") {
    ReadRec r;
    r.name = name; r.bc = bc; r.umi = umi;
    r.tid = 0; r.beg = beg; r.end = end; r.len = end - beg + 1;
    r.file = file; r.strand = 1;
    r.sjs = sjs;
    r.sigs.assign(sjs.size(), 1);
    v.push_back(std::move(r));
}

int main() {
    // ------------------------------------------------ interval algebra ----
    {
        IntervalVec a{{{10, 20}}}, b{{{15, 25}}};
        CHECK(intervalIntersection(a, b) == 6, "intersection [10,20] n [15,25] == 6");
        CHECK(intervalUnion(a, b) == 16, "union == 16");
        CHECK(chainDistance(a, b) == 10, "chain distance == 10");
        IntervalVec m{{{1, 5}, {4, 9}, {20, 22}}};
        CHECK(mergeIntervals(m).size() == 2, "mergeIntervals collapses overlap");
        CHECK(withinEditDistance("ACGTACGT", 8, "ACGTACGA", 8, 3), "UMI within 1 edit");
        CHECK(!withinEditDistance("AAAAAAAA", 8, "TTTTTTTT", 8, 3), "UMI far apart");
    }

    // ------------------------------------------------ GTF -----------------
    {
        std::ofstream f("/tmp/bc_test.gtf", std::ios::trunc);
        f << kGtf;
    }
    Annotation anno = loadGtf("/tmp/bc_test.gtf");
    const ChromAnno* ca = anno.find("chr1");
    CHECK(ca != nullptr, "chr1 present in annotation");
    CHECK(ca->tx_exons.size() == 3, "three transcripts parsed");
    CHECK(ca->tx_sjs.at("G1|T1").size() == 2, "T1 has two junctions");
    CHECK(ca->tx_sjs.at("G1|T1")[0] == (Interval{{201, 299}}), "T1 first junction");
    CHECK(ca->tx_sjs.count("G1|T2") == 1 && ca->tx_sjs.at("G1|T2").size() == 1,
          "T2 has one junction");
    CHECK(ca->tx_sjs.count("G2|T3") == 0, "T3 is single exon");
    CHECK(ca->gene_span.at("G1") == (Interval{{100, 600}}), "G1 span");

    // ------------------------------------------------ records -------------
    const IntervalVec chainT1{{{201, 299}, {401, 499}}};
    const IntervalVec chainT2{{{201, 499}}};
    const IntervalVec chainISM{{{201, 299}}};
    const IntervalVec chainNovel{{{201, 299}, {401, 520}}};

    std::vector<ReadRec> reads;
    for (int i = 0; i < 5; ++i) addRead(reads, "fsm1_" + std::to_string(i), 100, 600, chainT1);
    for (int i = 0; i < 3; ++i) addRead(reads, "fsm2_" + std::to_string(i), 100, 600, chainT2);
    for (int i = 0; i < 4; ++i) addRead(reads, "ism_" + std::to_string(i), 150, 400, chainISM);
    for (int i = 0; i < 6; ++i) addRead(reads, "nov_" + std::to_string(i), 100, 600, chainNovel);
    for (int i = 0; i < 2; ++i) addRead(reads, "se_" + std::to_string(i), 1010, 1090, {});

    {
        RecordWriter w("/tmp/bc_test_s0.rec");
        for (const auto& r : reads) w.append(r);
    }
    {   // round trip
        RecordReader rd("/tmp/bc_test_s0.rec");
        ReadRec r;
        int n = 0;
        bool ok = true;
        while (rd.next(r)) {
            if (r.name != reads[n].name || r.sjs != reads[n].sjs ||
                r.beg != reads[n].beg) ok = false;
            ++n;
        }
        CHECK(n == static_cast<int>(reads.size()) && ok, "record round trip");
    }

    // ------------------------------------------------ merge / grouping ----
    MergeResult merged = mergeAndGroup({{"/tmp/bc_test_s0.rec"}}, "/tmp/bc_test_merged.rec");
    CHECK(merged.n_records == static_cast<long long>(reads.size()), "all records merged");
    CHECK(merged.groups.size() == 2, "two clusters (chr1:100-600 and chr1:1010-1090)");
    CHECK(merged.groups[0].n_reads == 18 && merged.groups[1].n_reads == 2,
          "cluster sizes 18 / 2");

    Options opt;
    opt.mode = Mode::Bulk;
    opt.sj_support = 2;
    opt.graph_distance = 60;
    opt.min_count = 0;

    ColumnSpace cols;
    cols.init(Mode::Bulk, {"sample0"}, {});

    RecordReader reader("/tmp/bc_test_merged.rec");

    // ---- spliced cluster ----
    GroupData g0 = loadGroup(reader, merged.groups[0], "chr1", Mode::Bulk);
    CHECK(g0.read_sjs.size() == 18 && g0.single_exon.empty(),
          "18 spliced reads loaded");

    ClusterSet cs = buildClusters(g0);
    CHECK(cs.chain.size() == 4, "four distinct splice chains");

    std::vector<std::string> annoOrder;
    GroupAnno ga = sliceAnnotation(*ca, g0.span);
    for (const auto& kv : ga.me_sjs) annoOrder.push_back(kv.first);
    std::sort(annoOrder.begin(), annoOrder.end());
    CHECK(annoOrder.size() == 2, "two multi-exon models overlap the cluster");

    SpliceClasses sc;
    std::vector<TraceLine> trace;
    splitFsmIsmOthers(cs, g0, ga, buildAnnoIndex(ga, annoOrder), sc, trace);
    CHECK(sc.fsm["G1|T1"].size() == 5, "5 FSM reads for T1");
    CHECK(sc.fsm["G1|T2"].size() == 3, "3 FSM reads for T2");
    CHECK(sc.ism.size() == 1, "one ISM cluster");
    refineHighLow(cs, g0, opt.sj_support, sc);
    CHECK(sc.high.size() == 1, "novel chain reaches high confidence");

    GroupOutputs go = processGroup(g0, anno, opt, cols);
    CHECK(go.gtf.find("novel_isoform") != std::string::npos, "novel isoform emitted to GTF");
    CHECK(go.gtf.find("transcript_id \"T1\"") != std::string::npos, "T1 kept in GTF");
    CHECK(go.trace.find("\tFSM\t") != std::string::npos, "FSM lines in trace");
    CHECK(go.trace.find("\tISM\t") != std::string::npos, "ISM lines in trace");

    double txTotal = 0;
    for (const auto& r : go.tx_rows[0])
        for (const auto& v : r.vals) txTotal += v.second;
    CHECK(std::abs(txTotal - 18.0) < 1e-6,
          "transcript counts sum to 18 (got " + std::to_string(txTotal) + ")");

    double geneTotal = 0;
    for (const auto& r : go.gene_rows[0])
        for (const auto& v : r.vals) geneTotal += v.second;
    CHECK(std::abs(geneTotal - 18.0) < 1e-6,
          "gene counts sum to 18 (got " + std::to_string(geneTotal) + ")");

    // ---- single-exon cluster ----
    GroupData g1 = loadGroup(reader, merged.groups[1], "chr1", Mode::Bulk);
    CHECK(g1.single_exon.size() == 2, "two single-exon reads loaded");
    GroupOutputs go1 = processGroup(g1, anno, opt, cols);
    bool t3 = false;
    double t3n = 0;
    for (const auto& r : go1.tx_rows[0])
        if (r.name == "T3") { t3 = true; for (const auto& v : r.vals) t3n += v.second; }
    CHECK(t3 && std::abs(t3n - 2.0) < 1e-6, "single-exon reads land on T3");

    // ------------------------------------------------ EM sanity ----------
    {
        // one unambiguous read for A, four ambiguous reads compatible with A|B
        std::vector<std::vector<int>> rows{{0, 1}};
        std::vector<double> w{4.0}, prior{10.0, 0.0}, out;
        runEm(rows, w, prior, out);
        CHECK(out[0] + out[1] > 3.99 && out[0] + out[1] < 4.01,
              "EM conserves the ambiguous mass");
        CHECK(out[0] > out[1], "EM favours the transcript with prior evidence");
    }

    // ------------------------------------------------ table rewrite -------
    {
        TableSink sink("/tmp/bc_test_tx.sparse");
        std::vector<SparseRow> rows;
        SparseRow a; a.name = "TX1"; a.gene = "G1"; a.vals = {{0, 3.0}, {1, 1.0}};
        SparseRow b; b.name = "TX1"; b.gene = "G1"; b.vals = {{0, 2.0}};
        rows.push_back(a); rows.push_back(b);
        sink.write(rows);
        sink.close();

        RewriteOptions ro;
        rewriteTable("/tmp/bc_test_tx.sparse", "/tmp/bc_test_counts", {"s0", "s1"},
                     "transcript_id", ro);
        std::ifstream in("/tmp/bc_test_counts.txt");
        std::string header, line;
        std::getline(in, header);
        std::getline(in, line);
        CHECK(header == "transcript_id\tgene_id\ts0\ts1", "dense header");
        CHECK(line == "TX1\tG1\t5\t1", "duplicate rows summed (got '" + line + "')");
    }

    // ------------------------------------------------ read regions -------
    {
        // block reconstruction from span + junctions
        IntervalVec blk = readBlocks(Interval{{100, 600}}, chainT1);
        CHECK(blk.size() == 3 && blk[0] == (Interval{{100, 200}}) &&
              blk[1] == (Interval{{300, 400}}) && blk[2] == (Interval{{500, 600}}),
              "readBlocks rebuilds the three exons of T1");
        CHECK(readBlocks(Interval{{1010, 1090}}, {}).size() == 1,
              "readBlocks handles unspliced reads");

        GroupAnno ga2 = sliceAnnotation(*ca, Interval{{1, 5000}});
        CHECK(classifyRegion(blk, ga2, 0.5) == ReadRegion::Exonic, "spliced read is exonic");
        CHECK(classifyRegion({{{230, 280}}}, ga2, 0.5) == ReadRegion::Intronic,
              "read inside an intron of G1");
        CHECK(classifyRegion({{{3000, 3200}}}, ga2, 0.5) == ReadRegion::Intergenic,
              "read outside every gene");
        CHECK(classifyRegion({{{1020, 1080}}}, ga2, 0.5) == ReadRegion::Exonic,
              "read on the single-exon gene G2");
        CHECK(classifyRegion({{{590, 1010}}}, ga2, 0.5) == ReadRegion::Intergenic,
              "read straddling two loci but mostly between them");

        // overlapping genes must be reported as ambiguous
        {
            std::ofstream f("/tmp/bc_test_overlap.gtf", std::ios::trunc);
            f << "chr9\tsrc\texon\t100\t900\t.\t+\t.\tgene_id \"GA\"; transcript_id \"TA\";\n"
              << "chr9\tsrc\texon\t200\t800\t.\t-\t.\tgene_id \"GB\"; transcript_id \"TB\";\n";
        }
        Annotation ov = loadGtf("/tmp/bc_test_overlap.gtf");
        GroupAnno gov = sliceAnnotation(*ov.find("chr9"), Interval{{1, 2000}});
        CHECK(classifyRegion({{{300, 700}}}, gov, 0.5) == ReadRegion::Ambiguous,
              "read inside two overlapping genes is ambiguous");

        // the four categories must partition the reads of a cluster
        GroupOutputs gr = processGroup(g0, anno, opt, cols);
        double regTotal = 0;
        for (const auto& r : gr.region_rows[0])
            for (const auto& v : r.vals) regTotal += v.second;
        CHECK(std::abs(regTotal - 18.0) < 1e-6,
              "region counts partition the 18 reads (got " + std::to_string(regTotal) + ")");
        bool allExonic = false;
        for (const auto& r : gr.region_rows[0])
            if (r.name == std::string("Exonic") && r.vals[0].second == 18.0) allExonic = true;
        CHECK(allExonic, "all 18 reads of the G1 cluster are exonic");
    }

    // ------------------------------------------------ single-cell path ----
    {
        std::vector<ReadRec> screads;
        // cell1: 3 FSM reads for T1 with distinct UMIs, plus one 1-edit neighbour
        addRead(screads, "r1", 100, 600, chainT1, 0, "CELL1", "AAAAAAAA");
        addRead(screads, "r2", 100, 600, chainT1, 0, "CELL1", "CCCCCCCC");
        addRead(screads, "r3", 100, 600, chainT1, 0, "CELL1", "GGGGGGGG");
        addRead(screads, "r4", 100, 600, chainT1, 0, "CELL1", "GGGGGGGA");
        // cell2: 2 FSM reads for T2
        addRead(screads, "r5", 100, 600, chainT2, 0, "CELL2", "TTTTTTTT");
        addRead(screads, "r6", 100, 600, chainT2, 0, "CELL2", "TTTTTTAA");
        {
            RecordWriter w("/tmp/bc_test_sc.rec");
            for (const auto& r : screads) w.append(r);
        }
        MergeResult m2 = mergeAndGroup({{"/tmp/bc_test_sc.rec"}}, "/tmp/bc_test_sc_merged.rec");
        CHECK(m2.groups.size() == 1, "sc: one cluster");

        Options so;
        so.mode = Mode::SC;
        so.sj_support = 2;

        ColumnSpace sccols;
        sccols.init(Mode::SC, {"sample0"}, {{"CELL1", "CELL2"}});
        CHECK(sccols.tables() == 1 && sccols.columns(0) == 2, "sc: one table, two cells");
        CHECK(sccols.columnOf(0, "CELL2") == 1, "sc: barcode -> column index");

        RecordReader r2("/tmp/bc_test_sc_merged.rec");
        GroupData gs = loadGroup(r2, m2.groups[0], "chr1", Mode::SC);
        CHECK(gs.size() == 6, "sc: reads keyed by barcode-UMI");

        GroupOutputs so1 = processGroup(gs, anno, so, sccols);
        std::map<std::string, std::map<int, double>> byTx;
        for (const auto& row : so1.tx_rows[0])
            for (const auto& v : row.vals) byTx[row.name][v.first] += v.second;
        CHECK(byTx["T1"][0] == 4.0 && byTx["T1"].count(1) == 0,
              "sc: T1 counts land only in CELL1");
        CHECK(byTx["T2"][1] == 2.0, "sc: T2 counts land in CELL2");
        CHECK(so1.trace.find("CELL1") != std::string::npos, "sc: trace carries barcode");

        // now with fuzzy UMI collapsing enabled
        so.umi_dedup = true;
        so.umi_max_dist = 3;
        GroupOutputs so2 = processGroup(gs, anno, so, sccols);
        std::map<std::string, std::map<int, double>> byTx2;
        for (const auto& row : so2.tx_rows[0])
            for (const auto& v : row.vals) byTx2[row.name][v.first] += v.second;
        CHECK(byTx2["T1"][0] == 3.0, "sc: -k collapses the 1-edit UMI in CELL1");
        CHECK(byTx2["T2"][1] == 1.0, "sc: -k collapses the 2-edit UMI in CELL2");
    }

    std::cout << (g_fail ? "\nFAILURES: " : "\nALL TESTS PASSED (")
              << (g_fail ? std::to_string(g_fail) : std::string("0 failures)")) << "\n";
    return g_fail ? 1 : 0;
}
