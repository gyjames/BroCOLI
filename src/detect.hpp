// Splice-chain classification and novel-isoform detection.
//
// This is the algorithmic core and it is *identical* for bulk and single cell:
// the only single-cell specific step is the optional fuzzy UMI collapse, which
// is applied to the read lists before anything else looks at them.
#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common.hpp"
#include "gtf.hpp"
#include "group.hpp"

namespace brocoli {

// ------------------------------------------------------------- clustering --
struct ClusterSet {
    std::vector<IntervalVec> chain;                     // cluster -> splice chain
    std::vector<std::vector<std::string>> reads;        // cluster -> read keys
    std::vector<Interval> span;                         // cluster -> [min beg, max end]
};

inline ClusterSet buildClusters(const GroupData& g) {
    ClusterSet cs;
    std::map<IntervalVec, int> index;

    // Deterministic order: sort the read keys first so that cluster ids (and
    // therefore novel isoform names) do not depend on hash iteration order.
    std::vector<const std::string*> keys;
    keys.reserve(g.read_sjs.size());
    for (const auto& kv : g.read_sjs) keys.push_back(&kv.first);
    std::sort(keys.begin(), keys.end(),
              [](const std::string* a, const std::string* b) { return *a < *b; });

    for (const std::string* k : keys) {
        const IntervalVec& chain = g.read_sjs.at(*k);
        auto it = index.find(chain);
        int cid;
        if (it == index.end()) {
            cid = static_cast<int>(cs.chain.size());
            index.emplace(chain, cid);
            cs.chain.push_back(chain);
            cs.reads.emplace_back();
        } else {
            cid = it->second;
        }
        cs.reads[cid].push_back(*k);
    }

    cs.span.resize(cs.chain.size(), Interval{{0, 0}});
    for (size_t c = 0; c < cs.reads.size(); ++c) {
        int lo = INT_MAX, hi = INT_MIN;
        for (const auto& r : cs.reads[c]) {
            const Interval& iv = g.read_cov.at(r);
            lo = std::min(lo, iv[0]);
            hi = std::max(hi, iv[1]);
        }
        cs.span[c] = Interval{lo, hi};
    }
    return cs;
}

// ------------------------------------------------------------ UMI collapse --
// Greedy single-linkage collapse of reads whose UMIs are within `maxd` edits,
// applied per barcode. Replaces the edlib-based version of BroCOLI 1.x.
inline void collapseUmis(std::vector<std::string>& reads, const GroupData& g, int maxd) {
    if (reads.size() < 2) return;

    std::unordered_map<std::string, std::vector<size_t>> byBarcode;
    for (size_t i = 0; i < reads.size(); ++i) {
        auto it = g.read_bc.find(reads[i]);
        if (it != g.read_bc.end()) byBarcode[it->second].push_back(i);
    }

    std::vector<char> drop(reads.size(), 0);
    for (auto& kv : byBarcode) {
        auto& idx = kv.second;
        if (idx.size() < 2) continue;
        for (size_t a = 0; a < idx.size(); ++a) {
            if (drop[idx[a]]) continue;
            const std::string& ua = g.read_umi.at(reads[idx[a]]);
            for (size_t b = a + 1; b < idx.size(); ++b) {
                if (drop[idx[b]]) continue;
                const std::string& ub = g.read_umi.at(reads[idx[b]]);
                if (withinEditDistance(ua.data(), static_cast<int>(ua.size()),
                                       ub.data(), static_cast<int>(ub.size()), maxd)) {
                    drop[idx[b]] = 1;
                    ++g_stats.umi_collapsed;
                }
            }
        }
    }

    std::vector<std::string> kept;
    kept.reserve(reads.size());
    for (size_t i = 0; i < reads.size(); ++i)
        if (!drop[i]) kept.push_back(std::move(reads[i]));
    reads.swap(kept);
}

inline void collapseUmis(ClusterSet& cs, const GroupData& g, int maxd) {
    for (auto& reads : cs.reads) collapseUmis(reads, g, maxd);
}

// ---------------------------------------------------------- FSM / ISM / ? --
struct TraceLine {
    std::string read, category, isoform, gene;
};

struct SpliceClasses {
    std::unordered_map<std::string, std::vector<std::string>> fsm;  // tx key -> reads
    std::vector<int> ism;                                           // cluster ids
    std::vector<int> high;
    std::vector<int> low;
    std::unordered_map<int, std::string> high_strand;
};

// Is `part` a contiguous run inside `full`?
inline bool isContiguousSubchain(const IntervalVec& full, const IntervalVec& part) {
    if (part.empty() || part.size() > full.size()) return false;
    auto it = std::find(full.begin(), full.end(), part[0]);
    if (it == full.end()) return false;
    const size_t off = static_cast<size_t>(it - full.begin());
    if (off + part.size() > full.size()) return false;
    return std::equal(part.begin(), part.end(), full.begin() + off);
}

// Both FSM and ISM require the annotated chain to contain the cluster's first
// junction, so indexing annotations by junction turns the per-cluster scan over
// every transcript in the locus into a lookup. Lossless, and it is what makes
// dense loci (immunoglobulin, HLA, mitochondrial) tractable.
inline long long sjKey(const Interval& sj) {
    return (static_cast<long long>(sj[0]) << 32) | static_cast<uint32_t>(sj[1]);
}

struct AnnoIndex {
    std::vector<std::string> order;
    std::unordered_map<long long, std::vector<int>> by_first_sj;
};

inline AnnoIndex buildAnnoIndex(const GroupAnno& anno,
                                const std::vector<std::string>& order) {
    AnnoIndex ix;
    ix.order = order;
    ix.by_first_sj.reserve(order.size() * 4);
    for (size_t i = 0; i < order.size(); ++i) {
        auto it = anno.me_sjs.find(order[i]);
        if (it == anno.me_sjs.end()) continue;
        for (const auto& sj : it->second)
            ix.by_first_sj[sjKey(sj)].push_back(static_cast<int>(i));
    }
    return ix;
}

inline void splitFsmIsmOthers(const ClusterSet& cs, const GroupData& g,
                              const GroupAnno& anno, const AnnoIndex& ix,
                              SpliceClasses& out, std::vector<TraceLine>& trace) {
    std::vector<std::string> fsmHits, ismHits;

    for (size_t c = 0; c < cs.chain.size(); ++c) {
        const IntervalVec& chain = cs.chain[c];
        const auto& reads = cs.reads[c];
        if (reads.empty() || chain.empty()) continue;

        fsmHits.clear();
        ismHits.clear();

        auto bucket = ix.by_first_sj.find(sjKey(chain.front()));
        if (bucket != ix.by_first_sj.end())
        for (int ai : bucket->second) {
            const std::string& tx = ix.order[ai];
            const IntervalVec& ref = anno.me_sjs.at(tx);
            if (ref.size() == chain.size()) {
                if (ref == chain) fsmHits.push_back(tx);
            } else if (ref.size() > chain.size()) {
                if (isContiguousSubchain(ref, chain)) ismHits.push_back(tx);
            }
        }

        if (fsmHits.size() == 1) {
            auto& dst = out.fsm[fsmHits[0]];
            dst.insert(dst.end(), reads.begin(), reads.end());
            std::string gene, tx;
            splitTxKey(fsmHits[0], gene, tx);
            for (const auto& r : reads) trace.push_back({r, "FSM", tx, gene});

        } else if (!fsmHits.empty()) {
            // Several annotated transcripts share the chain: give each read to
            // the one whose annotated span is closest to the read's own span.
            for (const auto& r : reads) {
                const Interval& rc = g.read_cov.at(r);
                int best = 0, bestD = INT_MAX;
                for (size_t i = 0; i < fsmHits.size(); ++i) {
                    auto sp = anno.me_span.find(fsmHits[i]);
                    if (sp == anno.me_span.end()) continue;
                    const int d = std::abs(rc[0] - sp->second[0]) +
                                  std::abs(rc[1] - sp->second[1]);
                    if (d < bestD) { bestD = d; best = static_cast<int>(i); }
                }
                out.fsm[fsmHits[best]].push_back(r);
                std::string gene, tx;
                splitTxKey(fsmHits[best], gene, tx);
                trace.push_back({r, "FSM", tx, gene});
            }

        } else if (!ismHits.empty()) {
            out.ism.push_back(static_cast<int>(c));
        } else {
            out.high.push_back(static_cast<int>(c));   // provisional: refined below
        }
    }
}

// Splits the "unknown" clusters into high- and low-confidence sets. A cluster
// is high confidence when every one of its junctions is canonical in at least
// `support`+1 of its reads, the implied strand is consistent, and the chain is
// not a contiguous sub-chain of another high-confidence cluster.
inline void refineHighLow(const ClusterSet& cs, const GroupData& g,
                          int support, SpliceClasses& sc) {
    std::vector<int> others;
    others.swap(sc.high);

    std::vector<int> candidates;
    std::unordered_map<int, std::string> strand;
    std::vector<int> hits;

    for (int c : others) {
        const auto& reads = cs.reads[c];
        if (reads.size() <= 1) { sc.low.push_back(c); continue; }

        const size_t nsj = cs.chain[c].size();
        hits.assign(nsj, 0);
        int strandCode = 0;
        bool mixed = false;

        for (const auto& r : reads) {
            auto it = g.read_sigs.find(r);
            if (it == g.read_sigs.end()) continue;
            const auto& sig = it->second;
            for (size_t i = 0; i < sig.size() && i < nsj; ++i) {
                if (sig[i] == 1 || sig[i] == 2) {
                    ++hits[i];
                    if (strandCode && strandCode != sig[i]) mixed = true;
                    strandCode = sig[i];
                }
            }
        }

        bool allSupported = true;
        for (int h : hits) if (h <= support) { allSupported = false; break; }

        if (!allSupported || mixed || strandCode == 0) { sc.low.push_back(c); continue; }
        strand[c] = (strandCode == 1) ? "+" : "-";
        candidates.push_back(c);
    }

    // Demote clusters that are contained in a longer high-confidence chain.
    std::vector<char> demoted(cs.chain.size(), 0);
    for (int a : candidates)
        for (int b : candidates) {
            if (a == b || demoted[a]) continue;
            if (cs.chain[b].size() > cs.chain[a].size() &&
                isContiguousSubchain(cs.chain[b], cs.chain[a]))
                demoted[a] = 1;
        }

    for (int c : candidates) {
        if (demoted[c]) sc.low.push_back(c);
        else { sc.high.push_back(c); sc.high_strand[c] = strand[c]; }
    }
    std::sort(sc.high.begin(), sc.high.end());
    std::sort(sc.low.begin(), sc.low.end());
}

// Low-confidence clusters that are within `slack` of an annotated transcript
// with the same number of junctions are recycled as near-full-length matches.
inline void recycleLowConfidence(const ClusterSet& cs, const GroupAnno& anno,
                                 const std::vector<std::string>& annoOrder,
                                 SpliceClasses& sc, std::vector<TraceLine>& trace,
                                 int slack = 60) {
    if (sc.low.empty() || annoOrder.empty()) return;

    std::vector<int> remaining;
    std::vector<std::string> cand;
    std::vector<int> dist;

    for (int c : sc.low) {
        const IntervalVec& chain = cs.chain[c];
        cand.clear();
        dist.clear();
        for (const std::string& tx : annoOrder) {
            const IntervalVec& ref = anno.me_sjs.at(tx);
            if (ref.size() != chain.size()) continue;
            // Disjoint chains always score 0, so skipping them changes nothing.
            if (ref.front()[0] > chain.back()[1] || ref.back()[1] < chain.front()[0]) continue;
            const int d = chainDistance(ref, chain);
            if (d > 0 && d < slack) { cand.push_back(tx); dist.push_back(d); }
        }
        if (cand.empty()) { remaining.push_back(c); continue; }

        size_t pick = 0;
        if (cand.size() > 1) {
            // Prefer a transcript that already has full-length support,
            // otherwise the closest one.
            size_t bestFsm = 0;
            long bestCount = 0;
            for (size_t i = 0; i < cand.size(); ++i) {
                auto it = sc.fsm.find(cand[i]);
                const long n = (it == sc.fsm.end()) ? 0 : static_cast<long>(it->second.size());
                if (n > bestCount) { bestCount = n; bestFsm = i; }
            }
            if (bestCount > 0) pick = bestFsm;
            else pick = static_cast<size_t>(
                std::min_element(dist.begin(), dist.end()) - dist.begin());
        }

        auto& dst = sc.fsm[cand[pick]];
        dst.insert(dst.end(), cs.reads[c].begin(), cs.reads[c].end());
        std::string gene, tx;
        splitTxKey(cand[pick], gene, tx);
        for (const auto& r : cs.reads[c]) trace.push_back({r, "simFSM", tx, gene});
    }
    sc.low.swap(remaining);
}

// ------------------------------------------------- single-exon assignment --
struct SingleExonAssignment {
    std::unordered_map<std::string, std::vector<std::string>> gene_reads;
    std::unordered_map<std::string, std::vector<std::string>> tx_reads;
};

inline SingleExonAssignment assignSingleExon(const GroupData& g, const GroupAnno& anno,
                                             const ChromAnno& ca, const Options& opt,
                                             const SpliceClasses& sc,
                                             std::vector<TraceLine>& trace) {
    SingleExonAssignment out;
    if (anno.genes.empty() || g.single_exon.empty()) return out;

    std::vector<std::string> geneOrder;
    geneOrder.reserve(anno.genes.size());
    for (const auto& kv : anno.genes) geneOrder.push_back(kv.first);
    std::sort(geneOrder.begin(), geneOrder.end());

    std::vector<std::string> readOrder;
    readOrder.reserve(g.single_exon.size());
    for (const auto& kv : g.single_exon) readOrder.push_back(kv.first);
    std::sort(readOrder.begin(), readOrder.end());

    IntervalVec probe(1);
    std::vector<std::string> cand;
    std::vector<int> edge, overlap, nexon, tiebreak;

    // ---------------- reads -> gene ----------------
    for (const std::string& rk : readOrder) {
        const Interval& rex = g.single_exon.at(rk);
        probe[0] = rex;

        cand.clear(); edge.clear(); overlap.clear(); nexon.clear();
        for (const std::string& gene : geneOrder) {
            const Interval& gs = anno.genes.at(gene);
            if (rex[0] > gs[1] + opt.single_exon_edge ||
                rex[1] < gs[0] - opt.single_exon_edge) continue;

            auto ge = anno.gene_exons.find(gene);
            if (ge == anno.gene_exons.end()) continue;
            const int ov = intervalIntersection(ge->second, probe);
            if (ov <= 0) continue;

            int minDist = INT_MAX;
            auto g2t = ca.gene2tx.find(gene);
            if (g2t != ca.gene2tx.end())
                for (const std::string& tx : g2t->second) {
                    auto ex = ca.tx_exons.find(tx);
                    if (ex == ca.tx_exons.end()) continue;
                    minDist = std::min(minDist, intervalMinDistance(ex->second, rex));
                    if (minDist == 0) break;
                }

            cand.push_back(gene);
            edge.push_back(minDist);
            overlap.push_back(ov);
            nexon.push_back(static_cast<int>(ge->second.size()));
        }

        const int pick = pickCandidate(edge, overlap, nexon, {});
        if (pick >= 0) out.gene_reads[cand[pick]].push_back(rk);
    }

    if (opt.mode == Mode::SC && opt.umi_dedup)
        for (auto& kv : out.gene_reads) collapseUmis(kv.second, g, opt.umi_max_dist);

    // ---------------- reads -> transcript ----------------
    for (const auto& gkv : out.gene_reads) {
        auto g2t = ca.gene2tx.find(gkv.first);
        if (g2t == ca.gene2tx.end() || g2t->second.empty()) continue;
        const std::vector<std::string>& txs = g2t->second;

        for (const std::string& rk : gkv.second) {
            const Interval& rex = g.single_exon.at(rk);
            probe[0] = rex;

            if (txs.size() == 1) {
                out.tx_reads[txs[0]].push_back(rk);
                continue;
            }

            cand.clear(); edge.clear(); overlap.clear();
            nexon.clear(); tiebreak.clear();
            for (const std::string& tx : txs) {
                auto ex = ca.tx_exons.find(tx);
                if (ex == ca.tx_exons.end()) continue;
                const int ov = intervalIntersection(ex->second, probe);
                if (ov <= 0) continue;
                cand.push_back(tx);
                overlap.push_back(ov);
                nexon.push_back(static_cast<int>(ex->second.size()));
                edge.push_back(intervalMinDistance(ex->second, rex));
                auto f = sc.fsm.find(tx);
                tiebreak.push_back(f == sc.fsm.end() ? 0 : static_cast<int>(f->second.size()));
            }

            const int pick = pickCandidate(edge, overlap, nexon, tiebreak);
            if (pick >= 0) out.tx_reads[cand[pick]].push_back(rk);
        }
    }

    for (const auto& kv : out.tx_reads) {
        std::string gene, tx;
        splitTxKey(kv.first, gene, tx);
        for (const auto& r : kv.second) trace.push_back({r, "SE", tx, gene});
    }
    return out;
}

// ------------------------------------------------------- candidate graph ---
struct DistGraph {
    int n_anno = 0;
    std::vector<std::string> anno_name;       // node < n_anno
    std::vector<int> cluster_of;              // node >= n_anno -> cluster id (-1 otherwise)
    std::vector<int> count;                   // reads supporting the node
    std::vector<std::set<int>> adj;
    std::vector<std::string> novel_name;      // filled during output

    int nodes() const { return static_cast<int>(adj.size()); }
};

inline DistGraph buildGraph(const GroupAnno& anno,
                            const std::vector<std::string>& annoOrder,
                            const ClusterSet& cs, const SpliceClasses& sc,
                            int maxDist) {
    DistGraph gph;
    gph.n_anno = static_cast<int>(annoOrder.size());
    gph.anno_name = annoOrder;

    const int nNodes = gph.n_anno + static_cast<int>(sc.high.size());
    gph.adj.assign(nNodes, {});
    gph.count.assign(nNodes, 0);
    gph.cluster_of.assign(nNodes, -1);
    gph.novel_name.assign(nNodes, {});

    for (int i = 0; i < gph.n_anno; ++i) {
        auto it = sc.fsm.find(annoOrder[i]);
        gph.count[i] = (it == sc.fsm.end()) ? 0 : static_cast<int>(it->second.size());
    }
    for (size_t k = 0; k < sc.high.size(); ++k) {
        const int node = gph.n_anno + static_cast<int>(k);
        gph.cluster_of[node] = sc.high[k];
        gph.count[node] = static_cast<int>(cs.reads[sc.high[k]].size());
    }

    for (int i = 0; i < gph.n_anno; ++i) {
        const IntervalVec& a = anno.me_sjs.at(annoOrder[i]);
        for (size_t k = 0; k < sc.high.size(); ++k) {
            const IntervalVec& b = cs.chain[sc.high[k]];
            if (b.empty() || a.front()[0] > b.back()[1] || a.back()[1] < b.front()[0]) continue;
            const int node = gph.n_anno + static_cast<int>(k);
            const int d = chainDistance(a, b);
            if (d > 0 && d <= maxDist) {
                gph.adj[i].insert(node);
                gph.adj[node].insert(i);
            }
        }
    }
    for (size_t x = 0; x < sc.high.size(); ++x)
        for (size_t y = x + 1; y < sc.high.size(); ++y) {
            const IntervalVec& cx = cs.chain[sc.high[x]];
            const IntervalVec& cy = cs.chain[sc.high[y]];
            if (cx.empty() || cy.empty() ||
                cx.front()[0] > cy.back()[1] || cx.back()[1] < cy.front()[0]) continue;
            const int d = chainDistance(cx, cy);
            if (d > 0 && d <= maxDist) {
                const int nx = gph.n_anno + static_cast<int>(x);
                const int ny = gph.n_anno + static_cast<int>(y);
                gph.adj[nx].insert(ny);
                gph.adj[ny].insert(nx);
            }
        }
    return gph;
}

// Bron-Kerbosch with pivoting (iterative, as in BroCOLI 1.x).
inline std::vector<std::vector<int>> maximalCliques(const std::vector<std::set<int>>& adj,
                                                    size_t node_cap = 20000) {
    std::vector<std::vector<int>> out;
    const size_t n = adj.size();
    if (n == 0) return out;
    if (n > node_cap) {
        // Pathological locus: fall back to singletons rather than exploding.
        for (size_t i = 0; i < n; ++i) out.push_back({static_cast<int>(i)});
        return out;
    }

    std::set<int> cand, subg, ext;
    for (size_t i = 0; i < n; ++i) cand.insert(static_cast<int>(i));
    subg = cand;

    auto pivot = [&](const std::set<int>& s) {
        int best = s.empty() ? 0 : *s.begin();
        size_t bestDeg = 0;
        for (int v : s)
            if (adj[v].size() > bestDeg) { bestDeg = adj[v].size(); best = v; }
        return best;
    };

    std::vector<int> Q{-1};
    int p = pivot(subg);
    std::set_difference(cand.begin(), cand.end(), adj[p].begin(), adj[p].end(),
                        std::inserter(ext, ext.begin()));

    std::vector<std::tuple<std::set<int>, std::set<int>, std::set<int>>> stack;

    while (true) {
        if (!ext.empty()) {
            const int q = *ext.begin();
            ext.erase(ext.begin());
            cand.erase(q);
            Q.back() = q;

            std::set<int> subgQ;
            std::set_intersection(subg.begin(), subg.end(), adj[q].begin(), adj[q].end(),
                                  std::inserter(subgQ, subgQ.begin()));
            if (subgQ.empty()) {
                out.push_back(Q);
                continue;
            }
            std::set<int> candQ;
            std::set_intersection(cand.begin(), cand.end(), adj[q].begin(), adj[q].end(),
                                  std::inserter(candQ, candQ.begin()));
            if (candQ.empty()) continue;

            stack.emplace_back(subg, cand, ext);
            Q.push_back(-1);
            subg = std::move(subgQ);
            cand = std::move(candQ);
            p = pivot(subg);
            ext.clear();
            std::set_difference(cand.begin(), cand.end(), adj[p].begin(), adj[p].end(),
                                std::inserter(ext, ext.begin()));
        } else {
            if (stack.empty() || Q.empty()) break;
            Q.pop_back();
            std::tie(subg, cand, ext) = std::move(stack.back());
            stack.pop_back();
        }
    }
    return out;
}

struct DetectionResult {
    std::set<int> true_nodes;
    std::set<int> false_nodes;
};

inline bool cliqueLess(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.empty() || b.empty()) return b.empty() ? false : true;
    const int ma = *std::min_element(a.begin(), a.end());
    const int mb = *std::min_element(b.begin(), b.end());
    if (ma != mb) return ma < mb;
    return a.size() < b.size();
}

// Inside every maximal clique the annotated members are trusted; an unknown
// chain survives only if it is supported by at least as many reads as the
// weakest annotated member of the same clique.
inline DetectionResult detectTranscripts(std::vector<std::vector<int>>& cliques,
                                         const DistGraph& gph) {
    DetectionResult res;
    if (cliques.empty()) return res;

    std::set<int> known;
    for (int i = 0; i < gph.n_anno; ++i) known.insert(i);

    if (cliques.size() > 1) std::sort(cliques.begin(), cliques.end(), cliqueLess);

    for (const auto& cl : cliques) {
        std::set<int> cs(cl.begin(), cl.end());
        cs.erase(-1);

        std::set<int> knownHere, rest;
        std::set_intersection(cs.begin(), cs.end(), known.begin(), known.end(),
                              std::inserter(knownHere, knownHere.begin()));
        std::set_difference(cs.begin(), cs.end(), knownHere.begin(), knownHere.end(),
                            std::inserter(rest, rest.begin()));

        if (!knownHere.empty()) {
            int floorCount = INT_MAX;
            for (int n : knownHere) floorCount = std::min(floorCount, gph.count[n]);

            if (!rest.empty()) {
                res.true_nodes.insert(knownHere.begin(), knownHere.end());
                for (int n : rest) {
                    if (floorCount == 0 || gph.count[n] < floorCount) res.false_nodes.insert(n);
                    else res.true_nodes.insert(n);
                }
            } else {
                for (int n : knownHere)
                    if (gph.count[n] > 0) res.true_nodes.insert(n);
            }
        } else if (rest.size() == 1) {
            const int n = *rest.begin();
            if (gph.count[n] > 0) res.true_nodes.insert(n);
        } else {
            int best = 0;
            for (int n : rest) best = std::max(best, gph.count[n]);
            for (int n : rest) {
                if (gph.count[n] == best) res.true_nodes.insert(n);
                else res.false_nodes.insert(n);
            }
        }
    }
    for (int n : res.false_nodes) res.true_nodes.erase(n);
    return res;
}

}  // namespace brocoli
