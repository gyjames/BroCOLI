// One read cluster, start to finish. This function is the whole reason the two
// programs could be merged: bulk and single cell differ only in how a read maps
// to an output column, which is entirely hidden behind ColumnSpace.
#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.hpp"
#include "detect.hpp"
#include "gtf.hpp"
#include "group.hpp"
#include "output.hpp"
#include "quant.hpp"
#include "record.hpp"

namespace brocoli {

struct GroupOutputs {
    std::string gtf;
    std::string trace;
    std::vector<std::vector<SparseRow>> tx_rows;
    std::vector<std::vector<SparseRow>> gene_rows;
    std::vector<std::vector<SparseRow>> region_rows;
};

// Rebuild the aligned blocks of a read from its span and its splice junctions.
// A junction is stored as [last_exonic_base + 1, next_exonic_base - 1], so the
// blocks are exactly the complement of the junctions inside the span.
inline IntervalVec readBlocks(const Interval& span, const IntervalVec& sjs) {
    IntervalVec blocks;
    blocks.reserve(sjs.size() + 1);
    int cur = span[0];
    for (const auto& sj : sjs) {
        if (sj[0] - 1 >= cur) blocks.push_back(Interval{cur, sj[0] - 1});
        cur = sj[1] + 1;
    }
    if (span[1] >= cur) blocks.push_back(Interval{cur, span[1]});
    return blocks;
}

// Where does a read sit relative to the annotation?
//   Exonic     - at least `frac` of its aligned bases fall in the exons of
//                exactly one gene
//   Intronic   - inside exactly one gene locus, but not mostly in its exons
//   Ambiguous  - two or more genes make the same claim
//   Intergenic - no gene claims it
//
// This is an unstranded, purely positional call: it describes where the read
// landed on the genome and is independent of the FSM/ISM/novel assignment.
inline ReadRegion classifyRegion(const IntervalVec& blocks, const GroupAnno& anno,
                                 double frac) {
    int alen = 0;
    for (const auto& b : blocks) alen += b[1] - b[0] + 1;
    if (alen <= 0 || anno.genes.empty()) return ReadRegion::Intergenic;

    int nExonic = 0, nIntronic = 0;
    for (const auto& gk : anno.genes) {
        const IntervalVec locus{gk.second};
        const int lov = intervalIntersection(locus, blocks);
        if (lov <= 0) continue;

        auto ge = anno.gene_exons.find(gk.first);
        const int eov = (ge == anno.gene_exons.end())
                            ? 0
                            : intervalIntersection(ge->second, blocks);

        if (static_cast<double>(eov) / alen >= frac) ++nExonic;
        else if (static_cast<double>(lov) / alen >= frac) ++nIntronic;
        if (nExonic > 1) break;
    }

    if (nExonic > 1) return ReadRegion::Ambiguous;
    if (nExonic == 1) return ReadRegion::Exonic;
    if (nIntronic > 1) return ReadRegion::Ambiguous;
    if (nIntronic == 1) return ReadRegion::Intronic;
    return ReadRegion::Intergenic;
}

// Which annotated gene does a novel chain belong to? The gene whose pooled
// splice junctions overlap the chain most, provided the chain sits inside its
// locus; otherwise the closest enclosing gene; otherwise "NA".
inline std::string assignNovelGene(const IntervalVec& chain, const Interval& span,
                                   const GroupAnno& anno, const ChromAnno& ca) {
    std::string best;
    int bestOv = 0;

    for (const auto& gkv : anno.genes) {
        auto g2t = ca.gene2tx.find(gkv.first);
        if (g2t == ca.gene2tx.end()) continue;
        IntervalVec pooled;
        for (const std::string& tx : g2t->second) {
            auto sj = ca.tx_sjs.find(tx);
            if (sj != ca.tx_sjs.end())
                pooled.insert(pooled.end(), sj->second.begin(), sj->second.end());
        }
        if (pooled.empty()) continue;
        std::sort(pooled.begin(), pooled.end());
        pooled.erase(std::unique(pooled.begin(), pooled.end()), pooled.end());

        const int ov = intervalIntersection(pooled, chain);
        if (ov > bestOv || (ov == bestOv && ov > 0 && gkv.first < best)) {
            bestOv = ov;
            best = gkv.first;
        }
    }

    // A shared splice junction is decisive. BroCOLI 1.x additionally demanded
    // strict containment inside the gene locus, which sent every novel isoform
    // whose ends line up with the annotated gene to "NA" and silently removed
    // its reads from the gene-level counts.
    if (!best.empty()) return best;

    // No shared junction: fall back to genes that enclose the novel transcript.
    std::vector<std::string> enclosing;
    for (const auto& gkv : anno.genes)
        if (span[0] >= gkv.second[0] && span[1] <= gkv.second[1])
            enclosing.push_back(gkv.first);
    if (enclosing.empty()) return "NA";
    if (enclosing.size() == 1) return enclosing[0];

    std::sort(enclosing.begin(), enclosing.end());
    std::string pick = enclosing[0];
    int bestD = INT_MAX;
    for (const auto& gname : enclosing) {
        const Interval& gs = anno.genes.at(gname);
        const int d = std::min(std::abs(span[1] - gs[1]), std::abs(gs[0] - span[0]));
        if (d < bestD) { bestD = d; pick = gname; }
    }
    return pick;
}

inline GroupOutputs processGroup(const GroupData& g, const Annotation& annotation,
                                 const Options& opt, const ColumnSpace& cols) {
    GroupOutputs outp;
    outp.tx_rows.resize(cols.tables());
    outp.gene_rows.resize(cols.tables());
    outp.region_rows.resize(cols.tables());

    static const ChromAnno kEmptyChrom;
    const ChromAnno* caPtr = annotation.find(g.chrom);
    const ChromAnno& ca = caPtr ? *caPtr : kEmptyChrom;
    const GroupAnno anno = caPtr ? sliceAnnotation(ca, g.span) : GroupAnno{};

    std::vector<TraceLine> trace;
    trace.reserve(g.size());

    std::vector<std::string> annoOrder;
    annoOrder.reserve(anno.me_sjs.size());
    for (const auto& kv : anno.me_sjs) annoOrder.push_back(kv.first);
    std::sort(annoOrder.begin(), annoOrder.end());

    // ------------------------------------------------ spliced reads -------
    ClusterSet cs;
    SpliceClasses sc;
    if (!g.read_sjs.empty()) {
        cs = buildClusters(g);
        if (opt.mode == Mode::SC && opt.umi_dedup) collapseUmis(cs, g, opt.umi_max_dist);

        const AnnoIndex annoIndex = buildAnnoIndex(anno, annoOrder);
        splitFsmIsmOthers(cs, g, anno, annoIndex, sc, trace);
        refineHighLow(cs, g, opt.sj_support, sc);
        recycleLowConfidence(cs, anno, annoOrder, sc, trace);

        for (const auto& kv : sc.fsm) g_stats.fsm += static_cast<long long>(kv.second.size());
        for (int c : sc.ism)  g_stats.ism  += static_cast<long long>(cs.reads[c].size());
        for (int c : sc.high) g_stats.high += static_cast<long long>(cs.reads[c].size());
        for (int c : sc.low)  g_stats.low  += static_cast<long long>(cs.reads[c].size());
    }

    // ------------------------------------------------ single-exon reads ---
    const SingleExonAssignment se =
        g.single_exon.empty() ? SingleExonAssignment{}
                              : assignSingleExon(g, anno, ca, opt, sc, trace);

    // ------------------------------------------------ detection -----------
    const DistGraph gph = buildGraph(anno, annoOrder, cs, sc, opt.graph_distance);
    std::vector<std::vector<int>> cliques = maximalCliques(gph.adj);
    const DetectionResult det = detectTranscripts(cliques, gph);

    // ------------------------------------------------ build the model -----
    QuantModel qm;
    const bool skipNovel = (g.chrom == "chrM" || g.chrom == "MT" || g.chrom == "M");
    int novelCounter = 0;

    for (int node : det.true_nodes) {
        if (node < gph.n_anno) {
            const std::string& key = gph.anno_name[node];
            auto sj = anno.me_sjs.find(key);
            if (sj == anno.me_sjs.end()) continue;
            std::string gene, tx;
            splitTxKey(key, gene, tx);

            qm.node_to_tx[node] = static_cast<int>(qm.tx.size());
            qm.tx.push_back(key);
            qm.gene.push_back(gene);
            qm.chain.push_back(sj->second);

            auto ex = ca.tx_exons.find(key);
            auto st = ca.tx_strand.find(key);
            if (ex != ca.tx_exons.end())
                appendAnnotatedGtf(outp.gtf, g.chrom, key, ex->second,
                                   st == ca.tx_strand.end() ? "." : st->second);
        } else {
            if (skipNovel) continue;
            const int cid = gph.cluster_of[node];
            if (cid < 0) continue;
            ++novelCounter;
            const std::string name = g.chrom + "-novel-" + g.gid_str + "-" +
                                     std::to_string(node) + "-" + std::to_string(novelCounter);
            const IntervalVec& chain = cs.chain[cid];
            const Interval span = cs.span[cid];
            const std::string gene = anno.genes.empty()
                                         ? std::string("NA")
                                         : assignNovelGene(chain, span, anno, ca);

            std::string strand;
            if (gene != "NA") {
                auto gs = ca.gene_strand.find(gene);
                if (gs != ca.gene_strand.end()) strand = gs->second;
            }
            if (strand.empty()) {
                auto hs = sc.high_strand.find(cid);
                strand = (hs == sc.high_strand.end()) ? "." : hs->second;
            }

            qm.node_to_tx[node] = static_cast<int>(qm.tx.size());
            qm.tx.push_back(name);
            qm.gene.push_back(gene);
            qm.chain.push_back(chain);
            appendNovelGtf(outp.gtf, g.chrom, name, gene, span, chain, strand);

            for (const auto& r : cs.reads[cid]) trace.push_back({r, "novel", name, gene});
        }
    }

    addIsmRows(qm, cs, sc.ism);
    addFalseNodeRows(qm, gph, det, cs);

    // ------------------------------------------------ per column ----------
    const size_t T = qm.tx.size();
    // table -> tx index -> column -> count
    std::vector<std::vector<std::unordered_map<int, double>>> txCount(cols.tables());
    for (int t = 0; t < cols.tables(); ++t) txCount[t].resize(T);

    auto columnOf = [&](const std::string& read, int& table, int& col) {
        auto fit = g.read_file.find(read);
        if (fit == g.read_file.end()) return false;
        table = cols.tableOf(fit->second);
        if (opt.mode == Mode::Bulk) {
            col = cols.columnOf(fit->second, std::string());
        } else {
            auto bit = g.read_bc.find(read);
            if (bit == g.read_bc.end()) return false;
            col = cols.columnOf(fit->second, bit->second);
        }
        return table >= 0 && col >= 0;
    };

    if (T > 0) {
        // Unambiguous evidence (FSM reads, and reads of accepted novel chains)
        // bucketed by (table, column).
        std::vector<std::vector<std::unordered_map<int, double>>> prior(cols.tables());
        for (int t = 0; t < cols.tables(); ++t) prior[t].resize(T);

        std::unordered_map<std::string, int> txIndex;
        txIndex.reserve(T * 2);
        for (size_t j = 0; j < T; ++j) txIndex.emplace(qm.tx[j], static_cast<int>(j));

        for (const auto& kv : sc.fsm) {
            auto ti = txIndex.find(kv.first);
            if (ti == txIndex.end()) continue;
            for (const auto& r : kv.second) {
                int t, c;
                if (columnOf(r, t, c)) prior[t][ti->second][c] += 1.0;
            }
        }
        for (int node : det.true_nodes) {
            if (node < gph.n_anno) continue;
            auto nt = qm.node_to_tx.find(node);
            if (nt == qm.node_to_tx.end()) continue;
            const int cid = gph.cluster_of[node];
            if (cid < 0) continue;
            for (const auto& r : cs.reads[cid]) {
                int t, c;
                if (columnOf(r, t, c)) prior[t][nt->second][c] += 1.0;
            }
        }

        // Row weights: how many reads each indicator row contributes to a column.
        std::vector<std::vector<std::unordered_map<int, double>>> rowWeight(cols.tables());
        for (int t = 0; t < cols.tables(); ++t) rowWeight[t].resize(qm.rows.size());

        for (size_t i = 0; i < qm.rows.size(); ++i) {
            const int cid = (qm.row_kind[i] == 0) ? qm.row_src[i]
                                                  : gph.cluster_of[qm.row_src[i]];
            if (cid < 0) continue;

            // trace: which models the reads of this cluster are compatible with
            std::string txList, geneList;
            {
                std::set<std::string> tset, gset;
                for (int j : qm.rows[i]) { tset.insert(txPart(qm.tx[j])); gset.insert(qm.gene[j]); }
                for (const auto& s : tset) { if (!txList.empty()) txList += ','; txList += s; }
                for (const auto& s : gset) { if (!geneList.empty()) geneList += ','; geneList += s; }
            }
            const char* cat = (qm.row_kind[i] == 0) ? "ISM" : "simNovel";

            for (const auto& r : cs.reads[cid]) {
                int t, c;
                if (!columnOf(r, t, c)) continue;
                rowWeight[t][i][c] += 1.0;
                trace.push_back({r, cat, txList, geneList});
            }
        }

        std::vector<double> weights(qm.rows.size(), 0.0);
        std::vector<double> priorVec(T, 0.0), add(T, 0.0);

        for (int t = 0; t < cols.tables(); ++t) {
            // Which columns of this table are actually touched by this cluster?
            std::set<int> active;
            for (size_t j = 0; j < T; ++j)
                for (const auto& kv : prior[t][j]) active.insert(kv.first);
            for (size_t i = 0; i < qm.rows.size(); ++i)
                for (const auto& kv : rowWeight[t][i]) active.insert(kv.first);

            for (int c : active) {
                for (size_t j = 0; j < T; ++j) {
                    auto it = prior[t][j].find(c);
                    priorVec[j] = (it == prior[t][j].end()) ? 0.0 : it->second;
                }
                for (size_t i = 0; i < qm.rows.size(); ++i) {
                    auto it = rowWeight[t][i].find(c);
                    weights[i] = (it == rowWeight[t][i].end()) ? 0.0 : it->second;
                }
                runEm(qm.rows, weights, priorVec, add);
                for (size_t j = 0; j < T; ++j) {
                    const double v = priorVec[j] + add[j];
                    if (v > 0.0) txCount[t][j][c] += v;
                }
            }
        }
    }

    // ------------------------------------------------ gene level ----------
    // Gene counts come from the multi-exon transcripts plus the single-exon
    // reads assigned at gene level (matching BroCOLI 1.x, which deliberately
    // does not let single-exon reads count twice).
    std::vector<std::unordered_map<std::string, std::unordered_map<int, double>>>
        geneCount(cols.tables());

    for (int t = 0; t < cols.tables(); ++t)
        for (size_t j = 0; j < T; ++j) {
            if (qm.gene[j] == "NA") continue;
            auto& dst = geneCount[t][qm.gene[j]];
            for (const auto& kv : txCount[t][j]) dst[kv.first] += kv.second;
        }

    for (const auto& gkv : se.gene_reads)
        for (const auto& r : gkv.second) {
            int t, c;
            if (columnOf(r, t, c)) geneCount[t][gkv.first][c] += 1.0;
        }

    // ------------------------------------------------ single-exon counts --
    std::unordered_map<std::string, int> txIndex;
    for (size_t j = 0; j < T; ++j) txIndex.emplace(qm.tx[j], static_cast<int>(j));

    std::vector<std::unordered_map<std::string, std::unordered_map<int, double>>>
        seTxCount(cols.tables());

    for (const auto& kv : se.tx_reads) {
        if (kv.second.empty()) continue;
        const bool known = txIndex.count(kv.first) > 0;
        if (!known) {
            auto ex = ca.tx_exons.find(kv.first);
            auto st = ca.tx_strand.find(kv.first);
            if (ex != ca.tx_exons.end())
                appendAnnotatedGtf(outp.gtf, g.chrom, kv.first, ex->second,
                                   st == ca.tx_strand.end() ? "." : st->second);
        }
        for (const auto& r : kv.second) {
            int t, c;
            if (!columnOf(r, t, c)) continue;
            if (known) txCount[t][txIndex[kv.first]][c] += 1.0;
            else seTxCount[t][kv.first][c] += 1.0;
        }
    }

    // ------------------------------------------------ read regions -------
    // Every read in the cluster is classified exactly once, so the four
    // categories partition the reads and sum to the cluster size.
    std::vector<std::array<std::unordered_map<int, double>, kNumRegions>>
        regionCount(cols.tables());
    {
        long long local[kNumRegions] = {0, 0, 0, 0};
        static const IntervalVec kNoSj;
        for (const auto& kv : g.read_cov) {
            auto sj = g.read_sjs.find(kv.first);
            const IntervalVec blocks =
                readBlocks(kv.second, sj == g.read_sjs.end() ? kNoSj : sj->second);
            const int r = static_cast<int>(classifyRegion(blocks, anno, opt.region_frac));
            ++local[r];
            int t, c;
            if (columnOf(kv.first, t, c)) regionCount[t][r][c] += 1.0;
        }
        for (int r = 0; r < kNumRegions; ++r) g_stats.region[r] += local[r];
    }

    // ------------------------------------------------ emit rows -----------
    for (int t = 0; t < cols.tables(); ++t) {
        auto& rows = outp.tx_rows[t];
        for (size_t j = 0; j < T; ++j) {
            if (txCount[t][j].empty()) continue;
            SparseRow row;
            row.name = txPart(qm.tx[j]);
            row.gene = qm.gene[j];
            row.vals.assign(txCount[t][j].begin(), txCount[t][j].end());
            std::sort(row.vals.begin(), row.vals.end());
            rows.push_back(std::move(row));
        }
        for (const auto& kv : seTxCount[t]) {
            std::string gene, tx;
            splitTxKey(kv.first, gene, tx);
            SparseRow row;
            row.name = tx;
            row.gene = gene;
            row.vals.assign(kv.second.begin(), kv.second.end());
            std::sort(row.vals.begin(), row.vals.end());
            outp.tx_rows[t].push_back(std::move(row));
        }
        for (const auto& kv : geneCount[t]) {
            SparseRow row;
            row.name = kv.first;
            row.vals.assign(kv.second.begin(), kv.second.end());
            std::sort(row.vals.begin(), row.vals.end());
            outp.gene_rows[t].push_back(std::move(row));
        }
        for (int r = 0; r < kNumRegions; ++r) {
            if (regionCount[t][r].empty()) continue;
            SparseRow row;
            row.name = regionName(r);
            row.vals.assign(regionCount[t][r].begin(), regionCount[t][r].end());
            std::sort(row.vals.begin(), row.vals.end());
            outp.region_rows[t].push_back(std::move(row));
        }
    }

    // ------------------------------------------------ trace ---------------
    outp.trace.reserve(trace.size() * 64);
    for (const auto& tl : trace) {
        auto f = g.read_file.find(tl.read);
        if (f == g.read_file.end()) continue;
        outp.trace += tl.read;
        outp.trace += '\t';
        outp.trace += tl.category;
        outp.trace += '\t';
        outp.trace += tl.isoform;
        outp.trace += '\t';
        outp.trace += tl.gene;
        if (opt.mode == Mode::SC) {
            auto b = g.read_bc.find(tl.read);
            outp.trace += '\t';
            outp.trace += (b == g.read_bc.end() ? "." : b->second);
        }
        outp.trace += '\t';
        outp.trace += std::to_string(f->second);
        outp.trace += '\n';
    }
    return outp;
}

}  // namespace brocoli
