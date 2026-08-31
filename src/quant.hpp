// Quantification.
//
// Two changes relative to BroCOLI 1.x:
//
//  * Eigen is gone. The indicator matrix is 0/1 and extremely sparse, so it is
//    stored as one index list per row. The old code called
//    conservativeResize() once per row, i.e. it reallocated and copied the
//    whole dense matrix for every ISM cluster.
//
//  * The matrix is built once per read cluster and then reused for every
//    column. 1.x rebuilt it for every barcode, which is where most of the
//    single-cell runtime went.
#pragma once

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.hpp"
#include "detect.hpp"

namespace brocoli {

struct QuantModel {
    std::vector<std::string> tx;                 // matrix columns
    std::vector<std::string> gene;               // gene per tx
    std::vector<IntervalVec> chain;              // splice chain per tx
    std::unordered_map<int, int> node_to_tx;     // graph node -> index into tx

    std::vector<std::vector<int>> rows;          // indicator rows (tx indices)
    std::vector<int> row_kind;                   // 0 = ISM cluster, 1 = false node
    std::vector<int> row_src;                    // cluster id / graph node id

    bool empty() const { return tx.empty(); }
};

// Rows coming from ISM clusters: a cluster is compatible with every transcript
// that contains its chain as a contiguous run. When several transcripts match,
// only those with the highest chain-coverage ratio are kept (>= 0.5).
inline void addIsmRows(QuantModel& qm, const ClusterSet& cs,
                       const std::vector<int>& ismClusters) {
    std::vector<int> hit;
    for (int c : ismClusters) {
        const IntervalVec& chain = cs.chain[c];
        if (chain.empty()) continue;

        hit.clear();
        double maxRatio = 0.0;
        for (size_t j = 0; j < qm.chain.size(); ++j) {
            const IntervalVec& ref = qm.chain[j];
            if (ref.size() <= chain.size()) continue;
            if (chain.front()[0] < ref.front()[0] || chain.back()[1] > ref.back()[1]) continue;
            if (!isContiguousSubchain(ref, chain)) continue;
            hit.push_back(static_cast<int>(j));
            maxRatio = std::max(maxRatio,
                                static_cast<double>(chain.size()) / ref.size());
        }
        if (hit.empty()) continue;

        if (hit.size() > 1 && maxRatio >= 0.5) {
            std::vector<int> filtered;
            for (int j : hit) {
                const double r = static_cast<double>(chain.size()) / qm.chain[j].size();
                if (r >= maxRatio) filtered.push_back(j);
            }
            if (!filtered.empty()) hit.swap(filtered);
        }

        qm.rows.push_back(hit);
        qm.row_kind.push_back(0);
        qm.row_src.push_back(c);
    }
}

// Rows coming from rejected ("false") candidate chains: their reads are
// redistributed over the accepted neighbours in the candidate graph.
inline void addFalseNodeRows(QuantModel& qm, const DistGraph& gph,
                             const DetectionResult& det, const ClusterSet& cs) {
    for (int fn : det.false_nodes) {
        std::vector<int> row;

        for (int nb : gph.adj[fn]) {
            auto it = qm.node_to_tx.find(nb);
            if (it != qm.node_to_tx.end() && det.true_nodes.count(nb)) row.push_back(it->second);
        }

        if (row.empty()) {                       // second ring
            for (int nb : gph.adj[fn])
                for (int nn : gph.adj[nb]) {
                    auto it = qm.node_to_tx.find(nn);
                    if (it != qm.node_to_tx.end() && det.true_nodes.count(nn))
                        row.push_back(it->second);
                }
        }

        if (row.empty() && !qm.tx.empty()) {     // fall back to the closest model
            const int cid = gph.cluster_of[fn];
            if (cid >= 0) {
                const IntervalVec& q = cs.chain[cid];
                int best = 0, bestD = INT_MAX;
                for (size_t j = 0; j < qm.chain.size(); ++j) {
                    const int d = intervalUnion(qm.chain[j], q) -
                                  intervalIntersection(qm.chain[j], q);
                    if (d < bestD) { bestD = d; best = static_cast<int>(j); }
                }
                row.push_back(best);
            }
        }
        if (row.empty()) continue;

        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        qm.rows.push_back(std::move(row));
        qm.row_kind.push_back(1);
        qm.row_src.push_back(fn);
    }
}

// Standard EM over the compatibility matrix. `weights[i]` is how many reads
// row i carries in the current column, `prior[j]` the unambiguous (full-length)
// evidence already assigned to transcript j. `out[j]` receives the number of
// ambiguous reads apportioned to j.
inline void runEm(const std::vector<std::vector<int>>& rows,
                  const std::vector<double>& weights,
                  const std::vector<double>& prior,
                  std::vector<double>& out) {
    const size_t n = prior.size();
    out.assign(n, 0.0);
    if (rows.empty() || n == 0) return;

    double totalW = 0.0;
    for (double w : weights) totalW += w;
    if (totalW <= 0.0) return;

    std::vector<double> p(n, 1.0 / static_cast<double>(n));
    std::vector<double> next(n, 0.0);
    const int maxIter = (n > 500) ? 4 : 10;

    for (int iter = 0; iter < maxIter; ++iter) {
        std::fill(out.begin(), out.end(), 0.0);

        for (size_t i = 0; i < rows.size(); ++i) {
            const double w = weights[i];
            if (w <= 0.0 || rows[i].empty()) continue;
            double denom = 0.0;
            for (int j : rows[i]) denom += p[j];
            if (denom <= 0.0) continue;
            const double scale = w / denom;
            for (int j : rows[i]) out[j] += p[j] * scale;
        }

        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) { next[j] = out[j] + prior[j]; sum += next[j]; }
        if (sum <= 0.0) break;

        double delta = 0.0;
        for (size_t j = 0; j < n; ++j) {
            next[j] /= sum;
            delta += std::fabs(next[j] - p[j]);
        }
        p.swap(next);
        if (delta <= 5e-2) break;
    }
}

}  // namespace brocoli
