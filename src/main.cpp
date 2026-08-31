// BroCOLI 2.0 - unified bulk / single-cell long-read isoform detection and
// quantification, reading BAM through htslib.
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <streambuf>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "bam_stage1.hpp"
#include "common.hpp"
#include "detect.hpp"
#include "group.hpp"
#include "gtf.hpp"
#include "log.hpp"
#include "output.hpp"
#include "pipeline.hpp"
#include "quant.hpp"
#include "record.hpp"
#include "tpool.hpp"

using namespace brocoli;

namespace {

const char* kBanner = R"(
     ____               ____   ___  _      ___
    | __ )  _ __ ___   / ___| / _ \| |    |_ _|
    |  _ \ | '__/ _ \ | |    | | | | |     | |
    | |_) || | | (_) || |___ | |_| | |___  | |
    |____/ |_|  \___/  \____| \___/|_____||___|
)";

// Logging lives in log.hpp so preBroCOLI uses the identical format.
using brolog::fmt;
using brolog::logLine;
using brolog::LogRedirect;

struct option kLongOpts[] = {
    {"input",               required_argument, nullptr, 'i'},
    {"bam",                 required_argument, nullptr, 'i'},
    {"fasta",               required_argument, nullptr, 'f'},
    {"gtf",                 required_argument, nullptr, 'g'},
    {"output",              required_argument, nullptr, 'o'},
    {"mode",                required_argument, nullptr, 'M'},
    {"SJDistance",          required_argument, nullptr, 'j'},
    {"support",             required_argument, nullptr, 'n'},
    {"mapq",                required_argument, nullptr, 'q'},
    {"single_exon_boundary",required_argument, nullptr, 'e'},
    {"graph_distance",      required_argument, nullptr, 'd'},
    {"thread",              required_argument, nullptr, 't'},
    {"barcode",             required_argument, nullptr, 'B'},
    {"umi",                 required_argument, nullptr, 'U'},
    {"umi_dedup",           no_argument,       nullptr, 'k'},
    {"umi_dist",            required_argument, nullptr, 1001},
    {"output_min_count",    required_argument, nullptr, 'r'},
    {"window",              required_argument, nullptr, 1002},
    {"region_frac",         required_argument, nullptr, 1005},
    {"mtx",                 no_argument,       nullptr, 1003},
    {"keep_tmp",            no_argument,       nullptr, 1004},
    {"help",                no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

void printUsage(const char* prog) {
    std::cout
        << "Usage: " << prog << " -i <bam|list|dir> -f <genome.fa> -g <anno.gtf> -o <outdir>\n\n"
        << "Required:\n"
        << "  -i, --input                 coordinate-sorted+indexed BAM, a .txt/.tsv list of\n"
        << "                              BAMs (one per line), or a directory of BAMs\n"
        << "  -f, --fasta                 genome FASTA (chromosome names must match the BAM)\n"
        << "  -o, --output                output directory\n\n"
        << "Recommended:\n"
        << "  -g, --gtf                   reference annotation in GTF format\n"
        << "  -M, --mode                  bulk | sc                       (default: bulk)\n\n"
        << "Common:\n"
        << "  -t, --thread                worker threads                  (default: 8)\n"
        << "  -q, --mapq                  minimum mapping quality         (default: 0)\n"
        << "  -j, --SJDistance            minimum intron length           (default: 18)\n"
        << "  -n, --support               reads supporting each SJ of a novel isoform (2)\n"
        << "  -d, --graph_distance        candidate graph threshold       (default: 60)\n"
        << "  -e, --single_exon_boundary  slack for single-exon reads     (default: 300)\n"
        << "  -r, --output_min_count      drop counts below this value    (default: 0)\n\n"
        << "Single cell (--mode sc):\n"
        << "  -B, --barcode               cell barcode tag                (default: CB)\n"
        << "  -U, --umi                   UMI tag                         (default: UB)\n"
        << "  -k, --umi_dedup             collapse near-identical UMIs\n"
        << "      --umi_dist              edit distance for -k            (default: 3)\n"
        << "      --mtx                   write MatrixMarket instead of dense TSV\n\n"
        << "Advanced:\n"
        << "      --window                BAM interval size for stage 1   (default: 5000000)\n"
        << "      --region_frac           exonic/intronic call threshold  (default: 0.5)\n"
        << "      --keep_tmp              keep intermediate files\n"
        << "  -h, --help\n";
}

std::vector<std::string> discoverInputs(const std::string& path, const std::string& outdir) {
    std::vector<std::string> bams;
    struct stat st;
    if (stat(path.c_str(), &st) != 0) {
        std::cerr << "[BroCOLI] cannot stat " << path << ": " << std::strerror(errno) << "\n";
        return bams;
    }

    if (S_ISREG(st.st_mode)) {
        if (endsWith(path, ".txt") || endsWith(path, ".tsv")) {
            std::ifstream in(path);
            std::string line;
            while (std::getline(in, line)) {
                line = trimCopy(line);
                if (line.empty() || line[0] == '#') continue;
                bams.push_back(line);
            }
        } else if (endsWith(path, ".bam") || endsWith(path, ".cram")) {
            bams.push_back(path);
        } else {
            std::cerr << "[BroCOLI] expected .bam/.cram/.txt/.tsv, got " << path << "\n";
        }
    } else if (S_ISDIR(st.st_mode)) {
        DIR* d = opendir(path.c_str());
        if (!d) {
            std::cerr << "[BroCOLI] cannot open directory " << path << "\n";
            return bams;
        }
        while (struct dirent* e = readdir(d)) {
            std::string name = e->d_name;
            if (name.empty() || name[0] == '.') continue;
            if (endsWith(name, ".bam") || endsWith(name, ".cram"))
                bams.push_back(joinPath(path, name));
        }
        closedir(d);
        std::sort(bams.begin(), bams.end());
    }

    std::ofstream explain(joinPath(outdir, "file_explain.txt"), std::ios::trunc);
    explain << "File\tFile_Path\n";
    for (size_t i = 0; i < bams.size(); ++i) {
        std::cout << "[BroCOLI] sample " << i << ": " << bams[i] << "\n";
        explain << i << '\t' << bams[i] << '\n';
    }
    return bams;
}

// Union of every sample's @SQ lines, in the order of the first sample.
RefTable buildRefTable(const std::vector<std::string>& bams) {
    RefTable refs;
    for (const auto& b : bams) {
        samFile* fp = sam_open(b.c_str(), "r");
        if (!fp) { std::cerr << "[BroCOLI] cannot open " << b << "\n"; std::exit(EXIT_FAILURE); }
        sam_hdr_t* h = sam_hdr_read(fp);
        if (!h) { std::cerr << "[BroCOLI] cannot read header of " << b << "\n"; std::exit(EXIT_FAILURE); }
        for (int i = 0; i < sam_hdr_nref(h); ++i)
            refs.add(sam_hdr_tid2name(h, i), sam_hdr_tid2len(h, i));
        sam_hdr_destroy(h);
        hts_close(fp);
    }
    return refs;
}

void ensureIndex(const std::string& bam, int threads) {
    samFile* fp = sam_open(bam.c_str(), "r");
    if (!fp) return;
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());
    if (idx) { hts_idx_destroy(idx); hts_close(fp); return; }
    hts_close(fp);
    std::cout << "[BroCOLI] no index for " << bam << ", building one...\n";
    if (sam_index_build3(bam.c_str(), nullptr, 0, threads) < 0) {
        std::cerr << "[BroCOLI] failed to index " << bam
                  << " - is it coordinate sorted?\n";
        std::exit(EXIT_FAILURE);
    }
}

}  // namespace

int main(int argc, char** argv) {
    // No sync_with_stdio(false) here: cout and the C stdout must stay on the
    // same buffer, otherwise the banner and the progress lines land in a
    // redirected log in whatever order the two buffers happen to be flushed.
    std::cout << kBanner << "         BroCOLI  Version: 2.0.0\n" << std::flush;

    Options opt;
    int c, idx = 0;
    while ((c = getopt_long(argc, argv, "i:s:f:g:o:M:j:n:q:e:d:t:B:U:r:kh",
                            kLongOpts, &idx)) != -1) {
        switch (c) {
            case 'i': case 's': opt.input = optarg; break;
            case 'f': opt.fasta = optarg; break;
            case 'g': opt.gtf = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 'M': opt.mode = (std::string(optarg) == "sc") ? Mode::SC : Mode::Bulk; break;
            case 'j': opt.sj_distance = std::atoi(optarg); break;
            case 'n': opt.sj_support = std::atoi(optarg); break;
            case 'q': opt.mapq = std::atoi(optarg); break;
            case 'e': opt.single_exon_edge = std::atoi(optarg); break;
            case 'd': opt.graph_distance = std::atoi(optarg); break;
            case 't': opt.threads = std::atoi(optarg); break;
            case 'B': opt.bc_tag = optarg; break;
            case 'U': opt.umi_tag = optarg; break;
            case 'k': opt.umi_dedup = true; break;
            case 'r': opt.min_count = std::atoi(optarg); break;
            case 1001: opt.umi_max_dist = std::atoi(optarg); break;
            case 1002: opt.window = std::atoi(optarg); break;
            case 1003: opt.mtx = true; break;
            case 1004: opt.keep_tmp = true; break;
            case 1005: opt.region_frac = std::atof(optarg); break;
            case 'h': printUsage(argv[0]); return 0;
            default: printUsage(argv[0]); return EXIT_FAILURE;
        }
    }
    if (opt.input.empty() || opt.fasta.empty() || opt.outdir.empty()) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    if (opt.threads < 1) opt.threads = 1;
    if (opt.mode == Mode::Bulk && opt.mtx) opt.mtx = false;

    // From here on every line printed by any translation unit is timestamped
    // and flushed in order. Installed after option parsing so that the banner
    // and -h output stay plain.
    LogRedirect log_redirect;

    if (!makeDir(opt.outdir)) return EXIT_FAILURE;
    const std::string tmpdir = joinPath(opt.outdir, "tmp");
    if (!makeDir(tmpdir)) return EXIT_FAILURE;

    const std::vector<std::string> bams = discoverInputs(opt.input, opt.outdir);
    if (bams.empty()) {
        std::cerr << "[BroCOLI] no input BAM found\n";
        return EXIT_FAILURE;
    }
    if (bams.size() > 65535) {
        std::cerr << "[BroCOLI] too many samples\n";
        return EXIT_FAILURE;
    }

    std::cout << "*****\n"
              << "Mode:                 " << (opt.mode == Mode::SC ? "single cell" : "bulk") << '\n'
              << "Input:                " << opt.input << " (" << bams.size() << " sample(s))\n"
              << "FASTA:                " << opt.fasta << '\n'
              << "GTF:                  " << opt.gtf << '\n'
              << "Output:               " << opt.outdir << '\n'
              << "Threads:              " << opt.threads << '\n'
              << "MAPQ:                 " << opt.mapq << '\n'
              << "SJ distance:          " << opt.sj_distance << '\n'
              << "SJ support:           " << opt.sj_support << '\n'
              << "Graph distance:       " << opt.graph_distance << '\n'
              << "Single exon boundary: " << opt.single_exon_edge << '\n'
              << "Min output count:     " << opt.min_count << '\n';
    if (opt.mode == Mode::SC)
        std::cout << "Barcode tag:          " << opt.bc_tag << '\n'
                  << "UMI tag:              " << opt.umi_tag << '\n'
                  << "UMI dedup:            " << (opt.umi_dedup ? "on" : "off") << '\n';
    std::cout << "*****\n";

    for (const auto& b : bams) ensureIndex(b, opt.threads);

    const FastaRef fasta = loadFasta(opt.fasta);
    const RefTable refs = buildRefTable(bams);

    ThreadPool pool(opt.threads);
    auto t0 = std::chrono::steady_clock::now();

    // ------------------------------------------------------------ stage 1 --
    std::cout << "[BroCOLI] stage 1: scanning BAM files\n";
    std::vector<std::vector<std::string>> chunks(bams.size());
    std::vector<std::set<std::string>> barcodes(bams.size());
    for (size_t s = 0; s < bams.size(); ++s) {
        const std::string sdir = joinPath(tmpdir, "sample_" + std::to_string(s));
        if (!makeDir(sdir)) return EXIT_FAILURE;
        Stage1Result r = runStage1Sample(bams[s], static_cast<int>(s), opt, refs,
                                         fasta, sdir, pool);
        chunks[s] = std::move(r.chunk_files);
        barcodes[s] = std::move(r.barcodes);
        std::string line = "[BroCOLI]   " + bams[s] + ": " +
                           std::to_string(r.n_records) + " usable reads";
        if (opt.mode == Mode::SC)
            line += ", " + std::to_string(barcodes[s].size()) + " barcodes";
        logLine(line);
    }

    auto t1 = std::chrono::steady_clock::now();
    logLine(fmt("[BroCOLI] stage 1 took %.2f s",
                std::chrono::duration<double>(t1 - t0).count()));

    // ------------------------------------------------------------- merge ---
    const std::string mergedPath = joinPath(tmpdir, "merged.rec");
    MergeResult merged = mergeAndGroup(chunks, mergedPath);
    if (!opt.keep_tmp)
        for (const auto& v : chunks)
            for (const auto& f : v) std::remove(f.c_str());

    // --------------------------------------------------------- annotation --
    const Annotation annotation = loadGtf(opt.gtf);

    // ------------------------------------------------------------ outputs --
    std::vector<std::string> sampleNames;
    for (const auto& b : bams) sampleNames.push_back(baseNoExt(b));

    ColumnSpace cols;
    cols.init(opt.mode, sampleNames, barcodes);

    TextSink gtfSink(joinPath(opt.outdir, "updated_annotations.gtf"), "");
    TextSink traceSink(joinPath(opt.outdir, "compatible_isoform.tsv"),
                       opt.mode == Mode::SC
                           ? "read_id\tcategory\tisoform_id\tgene_id\tbarcode\tfile"
                           : "read_id\tcategory\tisoform_id\tgene_id\tfile");

    std::vector<std::unique_ptr<TableSink>> txSinks, geneSinks, regionSinks;
    for (int t = 0; t < cols.tables(); ++t) {
        txSinks.emplace_back(new TableSink(
            joinPath(tmpdir, "tx" + cols.tag(t) + ".sparse")));
        geneSinks.emplace_back(new TableSink(
            joinPath(tmpdir, "gene" + cols.tag(t) + ".sparse")));
        regionSinks.emplace_back(new TableSink(
            joinPath(tmpdir, "region" + cols.tag(t) + ".sparse")));
    }

    // ------------------------------------------------------ group workers --
    std::vector<size_t> order(merged.groups.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return merged.groups[a].n_reads > merged.groups[b].n_reads;
    });

    std::cout << "[BroCOLI] stage 3: identification and quantification\n";
    std::cout << "[BroCOLI] Display the results of clusters with more than 100,000 reads\n";
    t1 = std::chrono::steady_clock::now();
    {
        TaskQueue queue(pool);
        for (size_t k : order) {
            queue.submit([&, k] {
                thread_local std::unique_ptr<RecordReader> reader;
                if (!reader) reader.reset(new RecordReader(mergedPath));

                const GroupIndexEntry& e = merged.groups[k];
                const std::string& chrom = refs.name(e.tid);
                GroupData g = loadGroup(*reader, e, chrom, opt.mode);
                if (g.size() == 0) return;

                GroupOutputs go = processGroup(g, annotation, opt, cols);

                if (!go.gtf.empty()) gtfSink.write(go.gtf);
                if (!go.trace.empty()) traceSink.write(go.trace);
                for (int t = 0; t < cols.tables(); ++t) {
                    txSinks[t]->write(go.tx_rows[t]);
                    geneSinks[t]->write(go.gene_rows[t]);
                    regionSinks[t]->write(go.region_rows[t]);
                }
                if (e.n_reads > 100000)
                    logLine(fmt("[BroCOLI]   cluster %s:%d-%d, %d reads done",
                                chrom.c_str(), e.low, e.high, e.n_reads));
            });
        }
        queue.wait();
    }
    auto t2 = std::chrono::steady_clock::now();
    logLine(fmt("[BroCOLI] stage 3 took %.2f s",
                std::chrono::duration<double>(t2 - t1).count()));

    gtfSink.close();
    traceSink.close();
    for (auto& s : txSinks) s->close();
    for (auto& s : geneSinks) s->close();
    for (auto& s : regionSinks) s->close();

    // ------------------------------------------------------- final tables --
    std::cout << "[BroCOLI] writing count tables\n";
    RewriteOptions ro;
    ro.mtx = opt.mtx;
    ro.min_count = opt.min_count;

    for (int t = 0; t < cols.tables(); ++t) {
        ro.with_gene_column = true;
        rewriteTable(txSinks[t]->path(),
                     joinPath(opt.outdir, "counts_transcript" + cols.tag(t)),
                     cols.columnNames(t), "transcript_id", ro);
        ro.with_gene_column = false;
        rewriteTable(geneSinks[t]->path(),
                     joinPath(opt.outdir, "counts_gene" + cols.tag(t)),
                     cols.columnNames(t), "gene_id", ro);

        // Region counts are always a small dense table and are never filtered.
        RewriteOptions rr;
        rr.mtx = false;
        rr.min_count = 0;
        rr.with_gene_column = false;
        rewriteTable(regionSinks[t]->path(),
                     joinPath(opt.outdir, "read_regions" + cols.tag(t)),
                     cols.columnNames(t), "region", rr);
    }

    if (!opt.keep_tmp) {
        for (int t = 0; t < cols.tables(); ++t) {
            std::remove(txSinks[t]->path().c_str());
            std::remove(geneSinks[t]->path().c_str());
            std::remove(regionSinks[t]->path().c_str());
        }
        std::remove(mergedPath.c_str());
    }

    std::cout << "-----------------------------------------------------\n"
              << "Summary:\n"
              << "Alignments seen:       " << g_stats.reads_in.load() << '\n'
              << "Alignments used:       " << g_stats.reads_kept.load() << '\n'
              << "  dropped (flags):     " << g_stats.drop_flag.load() << '\n'
              << "  dropped (mapq):      " << g_stats.drop_mapq.load() << '\n'
              << "  dropped (cigar):     " << g_stats.drop_cigar.load() << '\n';
    if (opt.mode == Mode::SC)
        std::cout << "  dropped (no bc/umi): " << g_stats.drop_nobc.load() << '\n'
                  << "UMI collapsed:         " << g_stats.umi_collapsed.load() << '\n';
    std::cout << "Single exon reads:     " << g_stats.single_exon.load() << '\n'
              << "FSM reads:             " << g_stats.fsm.load() << '\n'
              << "ISM reads:             " << g_stats.ism.load() << '\n'
              << "High confidence reads: " << g_stats.high.load() << '\n'
              << "Low confidence reads:  " << g_stats.low.load() << '\n'
              << "-- genomic region --\n";
    {
        long long tot = 0;
        for (int r = 0; r < kNumRegions; ++r) tot += g_stats.region[r].load();
        for (int r = 0; r < kNumRegions; ++r) {
            const long long v = g_stats.region[r].load();
            logLine(fmt("%-22s %lld (%.2f%%)", regionName(r), v,
                        tot ? 100.0 * static_cast<double>(v) / static_cast<double>(tot) : 0.0));
        }
    }
    logLine(fmt("Total wall time:       %.2f s",
                std::chrono::duration<double>(t2 - t0).count()));
    std::cout << "-----------------------------------------------------\n" << std::flush;
    return 0;
}
