// BroCOLI 2.0 - common types and helpers
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace brocoli {

using Interval    = std::array<int, 2>;
using IntervalVec = std::vector<Interval>;

enum class Mode { Bulk, SC };

// ---------------------------------------------------------------- options --
struct Options {
    std::string input;              // .bam | list.txt/.tsv | directory
    std::string fasta;
    std::string gtf;
    std::string outdir;

    Mode mode            = Mode::Bulk;
    int  sj_distance     = 18;      // -j  minimum intron length
    int  sj_support      = 2;       // -n  reads supporting every SJ of a novel isoform
    int  graph_distance  = 60;      // -d  candidate graph edge threshold
    int  threads         = 8;       // -t
    int  min_count       = 0;       // -r  minimum count written out
    int  single_exon_edge= 300;     // -e  slack when attaching a single-exon read to a gene
    int  mapq            = 0;       // -q
    int  window          = 5000000; // --window  BAM interval granularity for stage 1

    std::string bc_tag   = "CB";    // -B  barcode tag (2 chars, "CB:Z:" also accepted)
    std::string umi_tag  = "UB";    // -U  umi tag
    bool umi_dedup       = false;   // -k  fuzzy (edit distance) UMI collapsing
    int  umi_max_dist    = 3;
    bool mtx             = false;   // --mtx  write MatrixMarket instead of dense TSV (sc)
    bool keep_tmp        = false;

    // Fraction of a read's aligned bases that must fall in a feature before the
    // read is called exonic / intronic for that gene.
    double region_frac   = 0.5;     // --region_frac
};

// ------------------------------------------------------- read region ------
enum class ReadRegion { Exonic = 0, Intronic = 1, Intergenic = 2, Ambiguous = 3 };
inline constexpr int kNumRegions = 4;

inline const char* regionName(int r) {
    switch (r) {
        case 0: return "Exonic";
        case 1: return "Intronic";
        case 2: return "Intergenic";
        default: return "Ambiguous";
    }
}

// ------------------------------------------------------------------ stats --
struct Stats {
    std::atomic<long long> reads_in{0};
    std::atomic<long long> reads_kept{0};
    std::atomic<long long> drop_flag{0};
    std::atomic<long long> drop_mapq{0};
    std::atomic<long long> drop_cigar{0};
    std::atomic<long long> drop_nobc{0};
    std::atomic<long long> single_exon{0};
    std::atomic<long long> fsm{0};
    std::atomic<long long> ism{0};
    std::atomic<long long> high{0};
    std::atomic<long long> low{0};
    std::atomic<long long> umi_collapsed{0};
    std::atomic<long long> region[4] = {};   // Exonic / Intronic / Intergenic / Ambiguous
};
inline Stats g_stats;

// -------------------------------------------------------- interval algebra --
// NOTE: the semantics below are inherited verbatim from BroCOLI 1.x so that
// every distance threshold keeps its original meaning.

inline int intervalIntersection(const IntervalVec& A, const IntervalVec& B) {
    int res = 0;
    size_t i = 0, j = 0;
    while (i < A.size() && j < B.size()) {
        const int s = std::max(A[i][0], B[j][0]);
        const int e = std::min(A[i][1], B[j][1]);
        if (s <= e) res += e - s + 1;
        if (A[i][1] < B[j][1]) ++i; else ++j;
    }
    return res;
}

inline int intervalUnion(const IntervalVec& A, const IntervalVec& B) {
    IntervalVec v;
    v.reserve(A.size() + B.size());
    v.insert(v.end(), A.begin(), A.end());
    v.insert(v.end(), B.begin(), B.end());
    std::sort(v.begin(), v.end());
    int res = 0;
    for (size_t i = 0; i < v.size();) {
        int l = v[i][0], r = v[i][1];
        while (i + 1 < v.size() && r >= v[i + 1][0]) { r = std::max(r, v[i + 1][1]); ++i; }
        res += r - l + 1;
        ++i;
    }
    return res;
}

inline int nonOverlappingLen(const IntervalVec& A, const IntervalVec& B) {
    int total = 0;
    for (const auto& a : A) {
        bool ov = false;
        for (const auto& b : B)
            if (!(a[1] < b[0] || a[0] > b[1])) { ov = true; break; }
        if (!ov) total += a[1] - a[0] + 1;
    }
    for (const auto& b : B) {
        bool ov = false;
        for (const auto& a : A)
            if (!(b[1] < a[0] || b[0] > a[1])) { ov = true; break; }
        if (!ov) total += b[1] - b[0] + 1;
    }
    return total;
}

// Distance between two splice chains: symmetric difference restricted to the
// region where the two chains actually overlap.
inline int chainDistance(const IntervalVec& A, const IntervalVec& B) {
    return intervalUnion(A, B) - intervalIntersection(A, B) - nonOverlappingLen(A, B);
}

inline int intervalMinDistance(const IntervalVec& A, const Interval& x) {
    int lmin = INT_MAX, rmin = INT_MAX;
    for (const auto& a : A) {
        lmin = std::min(lmin, std::abs(x[0] - a[0]));
        rmin = std::min(rmin, std::abs(x[1] - a[1]));
    }
    return std::min(lmin, rmin);
}

inline IntervalVec mergeIntervals(IntervalVec v) {
    if (v.empty()) return {};
    std::sort(v.begin(), v.end(),
              [](const Interval& a, const Interval& b) { return a[0] < b[0]; });
    IntervalVec out;
    out.reserve(v.size());
    out.push_back(v[0]);
    for (size_t i = 1; i < v.size(); ++i) {
        if (v[i][0] <= out.back()[1]) out.back()[1] = std::max(out.back()[1], v[i][1]);
        else out.push_back(v[i]);
    }
    return out;
}

// --------------------------------------------------------------- cascade ---
// Shared tie-breaking cascade used when a read could belong to several genes
// or several transcripts:
//   1. closest boundary, 2. largest overlap, 3. prefer single-exon models,
//   4. prefer the model that already collected the most full-length reads.
// Returns an index into the candidate arrays, or -1 when there is no candidate.
inline int pickCandidate(const std::vector<int>& edge,
                         const std::vector<int>& overlap,
                         const std::vector<int>& nexon,
                         const std::vector<int>& tiebreak) {
    const size_t n = edge.size();
    if (n == 0) return -1;
    if (n == 1) return 0;

    std::vector<int> s1;
    const int minEdge = *std::min_element(edge.begin(), edge.end());
    for (size_t i = 0; i < n; ++i)
        if (edge[i] - minEdge < 10) s1.push_back(static_cast<int>(i));
    if (s1.size() == 1) return s1[0];

    int maxOv = INT_MIN;
    for (int i : s1) maxOv = std::max(maxOv, overlap[i]);
    std::vector<int> s2;
    for (int i : s1)
        if (maxOv - overlap[i] < 10) s2.push_back(i);
    if (s2.size() == 1) return s2[0];

    std::vector<int> s3;
    if (!nexon.empty())
        for (int i : s2)
            if (nexon[i] == 1) s3.push_back(i);
    if (s3.size() == 1) return s3[0];

    const std::vector<int>& base = s3.empty() ? s2 : s3;
    if (tiebreak.empty()) return base[0];

    int best = base[0], bestVal = tiebreak[base[0]];
    for (int i : base)
        if (tiebreak[i] > bestVal) { bestVal = tiebreak[i]; best = i; }
    return best;
}

// ----------------------------------------------------------------- string --
inline bool endsWith(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() &&
           s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

inline std::string trimCopy(std::string x) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    x.erase(x.begin(), std::find_if(x.begin(), x.end(), notSpace));
    x.erase(std::find_if(x.rbegin(), x.rend(), notSpace).base(), x.end());
    return x;
}

inline std::string joinPath(const std::string& dir, const std::string& file) {
    if (dir.empty()) return file;
    return dir.back() == '/' ? dir + file : dir + "/" + file;
}

inline std::string baseNoExt(std::string p) {
    size_t pos = p.find_last_of("/\\");
    if (pos != std::string::npos) p = p.substr(pos + 1);
    pos = p.find_last_of('.');
    if (pos != std::string::npos && pos != 0) p = p.substr(0, pos);
    return p;
}

inline bool dirExists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR);
}

inline bool makeDir(const std::string& path) {
    if (dirExists(path)) return true;
    if (mkdir(path.c_str(), 0755) == 0) return true;
    std::fprintf(stderr, "[BroCOLI] cannot create directory %s: %s\n",
                 path.c_str(), std::strerror(errno));
    return false;
}

// "gene|transcript" -> gene / transcript
inline void splitTxKey(const std::string& key, std::string& gene, std::string& tx) {
    const size_t p = key.find('|');
    if (p == std::string::npos) { gene = key; tx = key; }
    else { gene = key.substr(0, p); tx = key.substr(p + 1); }
}

inline std::string txPart(const std::string& key) {
    const size_t p = key.find('|');
    return p == std::string::npos ? key : key.substr(p + 1);
}

// Banded Levenshtein, early-exits as soon as the distance provably exceeds
// `maxd`. Replaces the edlib dependency of BroCOLI 1.x; for 10-16 bp UMIs this
// is several times faster than a generic aligner.
inline bool withinEditDistance(const char* a, int la, const char* b, int lb, int maxd) {
    if (std::abs(la - lb) > maxd) return false;
    std::vector<int> prev(lb + 1), cur(lb + 1);
    for (int j = 0; j <= lb; ++j) prev[j] = j;
    for (int i = 1; i <= la; ++i) {
        cur[0] = i;
        const int lo = std::max(1, i - maxd);
        const int hi = std::min(lb, i + maxd);
        for (int j = 1; j < lo; ++j) cur[j] = maxd + 1;
        for (int j = lo; j <= hi; ++j) {
            const int sub = prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
            const int del = prev[j] + 1;
            const int ins = (j > lo ? cur[j - 1] : maxd + 1) + 1;
            cur[j] = std::min(sub, std::min(del, ins));
        }
        for (int j = hi + 1; j <= lb; ++j) cur[j] = maxd + 1;
        int rowMin = maxd + 1;
        for (int j = lo; j <= hi; ++j) rowMin = std::min(rowMin, cur[j]);
        if (rowMin > maxd) return false;
        prev.swap(cur);
    }
    return prev[lb] <= maxd;
}

}  // namespace brocoli
