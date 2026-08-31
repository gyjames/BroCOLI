// preBroCOLI 2.0 - per-read barcode and UMI assignment.
//
// This is a faithful port of the 1.x logic. P0 deliberately changes no
// decision the algorithm makes, so that its matched FASTQ can be diffed
// against the old binary's output. The only edits are ones that cannot alter
// a result:
//
//   * counters are atomics with readable names (they were racy bare ints);
//   * the edit-distance scratch buffer is thread_local instead of a fresh
//     heap allocation on every one of the millions of calls;
//   * an EdlibAlignResult leak in the no-whitelist branch is freed;
//   * two out-of-bounds reads (subpattern_ends walk, UMI padding on a read
//     truncated mid-UMI) are guarded.
//
// Everything marked FIXME(P1) is a real bug that DOES change results. They
// are left alone here on purpose and fixed in the next phase, one at a time,
// so each change can be A/B'd against real data.
#ifndef PREBROCOLI_DEMUX_HPP
#define PREBROCOLI_DEMUX_HPP

#include <chrono>
#include <numeric>

#include "common.hpp"
#include "whitelist.hpp"

#include "edlib.h"
#include "bindings/cpp/WFAligner.hpp"

namespace prebrocoli {

using Clock = std::chrono::steady_clock;
using microsec = std::chrono::microseconds;

// Banded-ish edit distance with early exit, exactly as in 1.x. The scratch
// matrix is reused per thread; 1.x allocated and zeroed a fresh vector on
// every call, and this function is on the innermost loop.
inline unsigned int edit_distance(const std::string& s1, const std::string& s2,
                                  unsigned int& end, int max_editd) {
    const std::size_t len1 = s1.size() + 1, len2 = s2.size() + 1;
    const char* s1_c = s1.c_str();
    const char* s2_c = s2.c_str();

    static thread_local std::vector<unsigned int> dist_holder;
    if (dist_holder.size() < len1 * len2) dist_holder.resize(len1 * len2);

    dist_holder[0] = 0;
    for (unsigned int j = 1; j < len2; ++j) dist_holder[j] = j;
    for (unsigned int i = 1; i < len1; ++i) dist_holder[i * len2] = 0;

    int best = static_cast<int>(len2);
    end = static_cast<unsigned int>(len1 - 1);

    for (unsigned int j = 1; j < len2; ++j) {
        bool any_below_threshold = false;
        for (unsigned int i = 1; i < len1; ++i) {
            const int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1;
            if (sub == 0) {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
            } else {
                dist_holder[i * len2 + j] = std::min({
                    dist_holder[(i - 1) * len2 + j] + 1,
                    dist_holder[i * len2 + (j - 1)] + 1,
                    dist_holder[(i - 1) * len2 + (j - 1)] + 1});
            }
            if (static_cast<int>(dist_holder[i * len2 + j]) <= max_editd)
                any_below_threshold = true;
            if (j == (len2 - 1) && static_cast<int>(dist_holder[i * len2 + j]) < best) {
                best = static_cast<int>(dist_holder[i * len2 + j]);
                end = i;
            }
        }
        if (!any_below_threshold) return 100;
    }
    return static_cast<unsigned int>(best);
}

// Full edit distance on two short buffers (<= 12 bp), stack DP table.
inline int edit_distance_ptr(const char* s1, int n1, const char* s2, int n2) {
    int dp[13][13];
    for (int i = 0; i <= n1; ++i) dp[i][0] = i;
    for (int j = 0; j <= n2; ++j) dp[0][j] = j;
    for (int i = 1; i <= n1; ++i)
        for (int j = 1; j <= n2; ++j) {
            const int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            dp[i][j] = std::min({dp[i - 1][j] + 1, dp[i][j - 1] + 1,
                                 dp[i - 1][j - 1] + cost});
        }
    return dp[n1][n2];
}

inline const EdlibEqualityPair* iupacEqualities(int& n) {
    static const EdlibEqualityPair kEq[32] = {
        {'R', 'A'}, {'R', 'G'}, {'K', 'G'}, {'K', 'T'},
        {'S', 'G'}, {'S', 'C'}, {'Y', 'C'}, {'Y', 'T'},
        {'M', 'A'}, {'M', 'C'}, {'W', 'A'}, {'W', 'T'},
        {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
        {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
        {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
        {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
        {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
        {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}};
    n = 32;
    return kEq;
}

// --------------------------------------------------------------------- UMI --
inline std::string get_umi(const std::string& seq, const Pattern& search_pattern,
                           const std::vector<int>& read_to_subpatterns,
                           int umi_index, UmiIndex& umi_map) {
    if (umi_index < 0) return "";
    const size_t umi_start = static_cast<size_t>(read_to_subpatterns[umi_index]);
    const size_t umi_length = search_pattern[umi_index].second.size();

    if (umi_start >= seq.size()) return std::string(umi_length, 'N');

    // 1.x computed the pad length as (pattern length - umi_length), which is
    // always zero, so a read truncated mid-UMI produced a short string that
    // was then indexed up to umi_length - an out-of-bounds read.
    std::string UMI_init = seq.substr(umi_start, umi_length);
    if (UMI_init.size() < umi_length)
        UMI_init.append(umi_length - UMI_init.size(), 'N');

    if (umi_map.valid_umis.empty()) return UMI_init;

    uint32_t exact_match_val = 0;
    bool can_encode = true;
    for (size_t i = 0; i < umi_length; ++i) {
        const unsigned char base = lut[static_cast<unsigned char>(UMI_init[i])];
        if (base == 0xFF) { can_encode = false; break; }
        exact_match_val = (exact_match_val << 2) | base;
    }

    if (can_encode &&
        std::binary_search(umi_map.valid_umis.begin(), umi_map.valid_umis.end(),
                           exact_match_val)) {
        bump(g_stats.umi_exact);
        return UMI_init;
    }

    static thread_local std::vector<int> votes;
    if (votes.size() < umi_map.valid_umis.size()) votes.resize(umi_map.valid_umis.size());
    std::memset(votes.data(), 0, umi_map.valid_umis.size() * sizeof(int));

    const int k = umi_map.k;
    const int len = static_cast<int>(UMI_init.size());
    if (len < k) return UMI_init;

    for (int p = 0; p <= len - k; ++p) {
        uint32_t kmer_val = 0;
        bool valid_kmer = true;
        for (int j = 0; j < k; ++j) {
            const unsigned char base = lut[static_cast<unsigned char>(UMI_init[p + j])];
            if (base == 0xFF) { valid_kmer = false; break; }
            kmer_val = (kmer_val << 2) | base;
        }
        if (!valid_kmer) continue;
        if (kmer_val < umi_map.kmer_lut.size())
            for (int cand : umi_map.kmer_lut[kmer_val]) votes[cand]++;
    }

    int best_id = -1, min_dist = 999;
    const int threshold_dist = 3;
    for (size_t i = 0; i < umi_map.valid_umis.size(); ++i) {
        if (votes[i] <= 0) continue;
        char buffer[16];
        const uint32_t val = umi_map.valid_umis[i];
        static const char dec[] = "ACGT";
        for (size_t j = 0; j < umi_length; ++j) {
            const int shift = 2 * static_cast<int>(umi_length - 1 - j);
            buffer[j] = dec[(val >> shift) & 3];
        }
        const int d = edit_distance_ptr(UMI_init.data(), static_cast<int>(UMI_init.size()),
                                        buffer, static_cast<int>(umi_length));
        if (d < min_dist) { min_dist = d; best_id = static_cast<int>(i); }
    }

    if (best_id != -1 && min_dist <= threshold_dist) {
        bump(g_stats.umi_corrected);
        return UmiIndex::decode32(umi_map.valid_umis[static_cast<size_t>(best_id)],
                                  static_cast<int>(umi_length));
    }

    bump(g_stats.umi_uncorrected);
    return UMI_init;
}

// Walks the edlib path and records where each sub-pattern starts in the read.
inline bool mapSubpatterns(const EdlibAlignResult& result, int flank_start,
                           const std::vector<unsigned long>& subpattern_lengths,
                           std::vector<int>& read_to_subpatterns) {
    std::vector<unsigned long> subpattern_ends(subpattern_lengths.size());
    std::partial_sum(subpattern_lengths.begin(), subpattern_lengths.end(),
                     subpattern_ends.begin());

    read_to_subpatterns.clear();
    read_to_subpatterns.reserve(subpattern_ends.size() + 1);
    read_to_subpatterns.push_back(flank_start);

    int i_read = flank_start;
    long i_pattern = 0;
    size_t i_subpattern = 0;

    for (int a = 0; a < result.alignmentLength; ++a) {
        const unsigned char value = result.alignment[a];
        if (value != 1) i_read++;
        if (value != 2) i_pattern++;
        // 1.x indexed subpattern_ends without checking i_subpattern first.
        if (i_subpattern < subpattern_ends.size() &&
            i_pattern >= static_cast<long>(subpattern_ends[i_subpattern])) {
            read_to_subpatterns.push_back(i_read);
            i_subpattern++;
        }
    }
    return read_to_subpatterns.size() >= subpattern_ends.size() + 1;
}

// ---------------------------------------------------------------- magicseq --
inline Barcode get_magicseq_barcode(std::string& seq, const BarcodeSet* known_bcX,
                                    const BarcodeSet* known_bcY, const BarcodeSet* known_bcZ,
                                    int flank_max_editd, const Pattern& search_pattern,
                                    const WhiteListMap* WHITELIST, BarcodeUMIindex& buindex) {
    Clock::time_point t1, t2;

    Barcode barcode;
    barcode.editd = -1000;
    barcode.flank_editd = 100;
    barcode.unambiguous = false;

    int n_eq = 0;
    const EdlibEqualityPair* eq = iupacEqualities(n_eq);
    EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, eq, n_eq};

    std::string search_string;
    search_string.reserve(105);
    std::vector<unsigned long> subpattern_lengths;
    subpattern_lengths.reserve(search_pattern.size());
    for (const auto& pair : search_pattern) {
        search_string += pair.second;
        subpattern_lengths.push_back(pair.second.length());
    }

    if (seq.length() < search_string.length()) {
        bump(g_stats.read_too_short);
        return barcode;
    }

    t1 = Clock::now();
    EdlibAlignResult result = edlibAlign(search_string.c_str(),
                                         static_cast<int>(search_string.length()),
                                         seq.c_str(), static_cast<int>(seq.length()),
                                         edlibConf);
    t2 = Clock::now();
    bump(g_stats.t_adapter_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (result.status != EDLIB_STATUS_OK || result.numLocations == 0) {
        bump(g_stats.adapter_not_found);
        edlibFreeAlignResult(result);
        return barcode;
    }

    barcode.flank_editd = result.editDistance;
    barcode.flank_start = result.startLocations[0];
    barcode.flank_end = result.endLocations[0];

    std::vector<int> read_to_subpatterns;
    const bool mapped = mapSubpatterns(result, barcode.flank_start,
                                       subpattern_lengths, read_to_subpatterns);
    edlibFreeAlignResult(result);
    if (!mapped) {
        bump(g_stats.adapter_not_found);
        return barcode;
    }

    const int bcx_index = 1, bcy_index = 3, bcz_index = 5, umi_index = 6;
    const int edit_XY = 10;
    const int edit_Z = edit_XY - 2;

    std::string BarcodeX, BarcodeY, BarcodeZ;
    std::vector<std::string> BarcodeVecStrX, BarcodeVecStrY, BarcodeVecStrZ;

    // ------------------------------------------------------------ barcode X --
    t1 = Clock::now();
    std::string exact_bcx = seq.substr(
        static_cast<size_t>(read_to_subpatterns[bcx_index]),
        static_cast<size_t>(read_to_subpatterns[bcx_index + 1] - read_to_subpatterns[bcx_index]));
    int barcodeX_editd = 100;
    if (known_bcX->empty() || known_bcX->find(exact_bcx) != known_bcX->end()) {
        barcodeX_editd = 0;
        BarcodeX = exact_bcx;
        BarcodeVecStrX.push_back(BarcodeX);
    } else {
        std::unordered_map<std::string, int> candidate_votes;
        std::string bcx_ext;
        bcx_ext.reserve(exact_bcx.size() + 2);
        bcx_ext.append("T"); bcx_ext += exact_bcx; bcx_ext.append("C");

        std::string kmer(magicseqXKMER, '\0');
        for (size_t p = 0; p + magicseqXKMER <= bcx_ext.size(); ++p) {
            std::memcpy(&kmer[0], bcx_ext.data() + p, magicseqXKMER);
            auto it = buindex.barcodeX_kmer_index.find(kmer);
            if (it != buindex.barcodeX_kmer_index.end())
                for (const std::string& id : it->second) candidate_votes[id]++;
        }
        // FIXME(P1): unordered_map iteration order + unstable sort means ties
        // are broken by hash order, so top-8 is not reproducible across builds.
        std::vector<std::pair<std::string, int>> cands(candidate_votes.begin(),
                                                       candidate_votes.end());
        std::sort(cands.begin(), cands.end(),
                  [](const std::pair<std::string, int>& a,
                     const std::pair<std::string, int>& b) { return a.second > b.second; });
        if (cands.size() > 8) cands.resize(8);

        int left_bound = std::max(read_to_subpatterns[bcx_index - 1], 0);
        int max_length = read_to_subpatterns[bcx_index + 2] - read_to_subpatterns[bcx_index - 1];
        std::string barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                             static_cast<size_t>(max_length));

        unsigned int editDistance, endDistance;
        for (const auto& cand : cands) {
            search_string.clear();
            search_string.append(search_pattern[bcx_index - 1].second);
            search_string.append(cand.first);
            search_string.append(search_pattern[bcx_index + 1].second);
            editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_XY);
            if (static_cast<int>(editDistance) == barcodeX_editd) {
                if (static_cast<int>(editDistance) < edit_XY)
                    BarcodeVecStrX.push_back(cand.first);
            } else if (static_cast<int>(editDistance) < barcodeX_editd &&
                       static_cast<int>(editDistance) < edit_XY) {
                BarcodeVecStrX.clear();
                BarcodeVecStrX.push_back(cand.first);
                barcodeX_editd = static_cast<int>(editDistance);
                BarcodeX = cand.first;
            }
        }

        if (BarcodeVecStrX.empty()) {
            left_bound = std::max(read_to_subpatterns[bcx_index], 0);
            max_length = read_to_subpatterns[bcx_index + 1] - read_to_subpatterns[bcx_index];
            barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                     static_cast<size_t>(max_length));
            for (const auto& cand : cands) {
                search_string.clear();
                search_string.append(cand.first);
                editDistance = edit_distance(barcode_seq, search_string, endDistance, 2);
                if (editDistance < 3) {
                    BarcodeVecStrX.push_back(cand.first);
                    barcodeX_editd = static_cast<int>(editDistance);
                    BarcodeX = cand.first;
                }
            }
        }
    }
    t2 = Clock::now();
    bump(g_stats.t_bc1_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (BarcodeVecStrX.empty() && BarcodeX.empty()) {
        bump(g_stats.bc1_unresolved);
        return barcode;
    }
    if (BarcodeVecStrX.size() == 1) exact_bcx = BarcodeVecStrX[0];

    // ------------------------------------------------------------ barcode Y --
    t1 = Clock::now();
    std::string exact_bcy = seq.substr(
        static_cast<size_t>(read_to_subpatterns[bcy_index]),
        static_cast<size_t>(read_to_subpatterns[bcy_index + 1] - read_to_subpatterns[bcy_index]));
    int barcodeY_editd = 100;
    if (known_bcY->empty() || known_bcY->find(exact_bcy) != known_bcY->end()) {
        barcodeY_editd = 0;
        BarcodeY = exact_bcy;
        BarcodeVecStrY.push_back(BarcodeY);
    } else {
        std::unordered_map<std::string, int> candidate_votes;
        std::string bcy_ext;
        bcy_ext.reserve(exact_bcy.size() + 2);
        bcy_ext.append("A"); bcy_ext += exact_bcy; bcy_ext.append("T");

        std::string kmer(magicseqXKMER, '\0');
        for (size_t p = 0; p + magicseqXKMER <= bcy_ext.size(); ++p) {
            std::memcpy(&kmer[0], bcy_ext.data() + p, magicseqXKMER);
            auto it = buindex.barcodeY_kmer_index.find(kmer);
            if (it != buindex.barcodeY_kmer_index.end())
                for (const std::string& id : it->second) candidate_votes[id]++;
        }
        std::vector<std::pair<std::string, int>> cands(candidate_votes.begin(),
                                                       candidate_votes.end());
        std::sort(cands.begin(), cands.end(),
                  [](const std::pair<std::string, int>& a,
                     const std::pair<std::string, int>& b) { return a.second > b.second; });
        if (cands.size() > 8) cands.resize(8);

        int left_bound = std::max(read_to_subpatterns[bcy_index - 1], 0);
        int max_length = read_to_subpatterns[bcy_index + 2] - read_to_subpatterns[bcy_index - 1];
        std::string barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                             static_cast<size_t>(max_length));

        unsigned int editDistance, endDistance;
        for (const auto& cand : cands) {
            search_string.clear();
            search_string.append(search_pattern[bcy_index - 1].second);
            search_string.append(cand.first);
            search_string.append(search_pattern[bcy_index + 1].second);
            editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_XY);
            if (static_cast<int>(editDistance) == barcodeY_editd) {
                if (static_cast<int>(editDistance) < edit_XY)
                    BarcodeVecStrY.push_back(cand.first);
            } else if (static_cast<int>(editDistance) < barcodeY_editd &&
                       static_cast<int>(editDistance) < edit_XY) {
                BarcodeVecStrY.clear();
                BarcodeVecStrY.push_back(cand.first);
                barcodeY_editd = static_cast<int>(editDistance);
                BarcodeY = cand.first;
            }
        }

        if (BarcodeVecStrY.empty()) {
            left_bound = std::max(read_to_subpatterns[bcy_index], 0);
            max_length = read_to_subpatterns[bcy_index + 1] - read_to_subpatterns[bcy_index];
            barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                     static_cast<size_t>(max_length));
            for (const auto& cand : cands) {
                search_string.clear();
                search_string.append(cand.first);
                editDistance = edit_distance(barcode_seq, search_string, endDistance, 2);
                if (editDistance < 3) {
                    BarcodeVecStrY.push_back(cand.first);
                    barcodeY_editd = static_cast<int>(editDistance);
                    BarcodeY = cand.first;
                }
            }
        }
    }
    t2 = Clock::now();
    bump(g_stats.t_bc2_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (BarcodeVecStrY.empty() && BarcodeY.empty()) {
        bump(g_stats.bc2_unresolved);
        return barcode;
    }
    if (BarcodeVecStrY.size() == 1) exact_bcy = BarcodeVecStrY[0];

    // ------------------------------------------------------------ barcode Z --
    // FIXME(P1): this branch scans the whole barcodeZ set linearly while
    // buindex.barcodeZ_kmer_index sits unused, and the rescue loop below is
    // dead code (the iterator is already at end()).
    t1 = Clock::now();
    std::string exact_bcz = seq.substr(
        static_cast<size_t>(read_to_subpatterns[bcz_index]),
        static_cast<size_t>(read_to_subpatterns[bcz_index + 1] - read_to_subpatterns[bcz_index]));
    int barcodeZ_editd = 100;
    if (known_bcZ->empty() || known_bcZ->find(exact_bcz) != known_bcZ->end()) {
        barcodeZ_editd = 0;
        BarcodeZ = exact_bcz;
        BarcodeVecStrZ.push_back(BarcodeZ);
    } else {
        const int left_bound = std::max(read_to_subpatterns[bcz_index - 1], 0);
        const int max_length = read_to_subpatterns[bcz_index + 1] - read_to_subpatterns[bcz_index - 1];
        const std::string barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                                   static_cast<size_t>(max_length));
        unsigned int editDistance, endDistance;
        auto itr = known_bcZ->begin();
        for (; itr != known_bcZ->end(); ++itr) {
            search_string.clear();
            search_string.append(search_pattern[bcz_index - 1].second);
            search_string.append(*itr);
            editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_Z);
            if (static_cast<int>(editDistance) == barcodeZ_editd) {
                if (static_cast<int>(editDistance) < edit_Z)
                    BarcodeVecStrZ.push_back(*itr);
            } else if (static_cast<int>(editDistance) < barcodeZ_editd &&
                       static_cast<int>(editDistance) < edit_Z) {
                BarcodeVecStrZ.clear();
                BarcodeVecStrZ.push_back(*itr);
                barcodeZ_editd = static_cast<int>(editDistance);
                BarcodeZ = *itr;
            }
        }
        // FIXME(P1): itr == end() here, so this rescue never runs.
        if (BarcodeVecStrZ.empty()) {
            for (; itr != known_bcZ->end(); ++itr) { /* unreachable in 1.x */ }
        }
    }
    t2 = Clock::now();
    bump(g_stats.t_bc3_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (BarcodeVecStrZ.empty() && BarcodeZ.empty()) {
        bump(g_stats.bc3_unresolved);
        return barcode;
    }
    if (BarcodeVecStrZ.size() == 1) exact_bcz = BarcodeVecStrZ[0];

    // ----------------------------------------------------- whitelist lookup --
    std::string Final_Barcode;
    Final_Barcode.reserve(24);
    t1 = Clock::now();
    if (!WHITELIST->empty()) {
        std::string whiteListString;
        whiteListString.reserve(24);
        std::unordered_set<std::string> Final24;
        std::string kmer(magicseqALLKMER, '\0');
        for (const auto& x : BarcodeVecStrX)
            for (const auto& y : BarcodeVecStrY)
                for (const auto& z : BarcodeVecStrZ) {
                    // FIXME(P1): clear() belongs outside the three loops; as
                    // written only the last (x,y,z) combination survives.
                    Final24.clear();
                    whiteListString.clear();
                    whiteListString.append(x);
                    whiteListString.append(y);
                    whiteListString.append(z);
                    for (size_t p = 0; p + magicseqALLKMER <= whiteListString.size(); ++p) {
                        std::memcpy(&kmer[0], whiteListString.data() + p, magicseqALLKMER);
                        auto it = buindex.whiteList_kmer_index.find(kmer);
                        if (it != buindex.whiteList_kmer_index.end())
                            for (const std::string& id : it->second) Final24.insert(id);
                    }
                }

        wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment,
                                        wfa::WFAligner::MemoryHigh);
        const int left_bound = std::max(read_to_subpatterns[bcx_index - 1], 0);
        const int max_length = read_to_subpatterns[bcz_index + 1] - read_to_subpatterns[bcx_index - 1];
        const std::string barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                                   static_cast<size_t>(max_length));
        for (const auto& fin : Final24) {
            search_string.clear();
            search_string.append(search_pattern[0].second); search_string.append(fin.substr(0, 8));
            search_string.append(search_pattern[2].second); search_string.append(fin.substr(8, 8));
            search_string.append(search_pattern[4].second); search_string.append(fin.substr(16, 8));
            aligner.alignEnd2End(search_string, barcode_seq);
            if (aligner.getAlignmentScore() > barcode.editd) {
                Final_Barcode = fin;
                barcode.editd = aligner.getAlignmentScore();
            }
        }
    } else {
        for (const auto& x : BarcodeVecStrX)
            for (const auto& y : BarcodeVecStrY)
                for (const auto& z : BarcodeVecStrZ) {
                    search_string.clear();
                    search_string.append(search_pattern[0].second); search_string.append(x);
                    search_string.append(search_pattern[2].second); search_string.append(y);
                    search_string.append(search_pattern[4].second); search_string.append(z);
                    search_string.append("NNNNNNNNNNNN");
                    search_string.append("TTTTTTTTTT");
                    EdlibAlignResult r2 = edlibAlign(search_string.c_str(),
                                                     static_cast<int>(search_string.length()),
                                                     seq.c_str(), static_cast<int>(seq.length()),
                                                     edlibConf);
                    if (r2.editDistance < barcode.editd) {
                        BarcodeX = x; BarcodeY = y; BarcodeZ = z;
                        barcode.editd = r2.editDistance;
                    }
                    edlibFreeAlignResult(r2);   // 1.x leaked this on every call
                }
    }
    t2 = Clock::now();
    bump(g_stats.t_whitelist_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (!Final_Barcode.empty()) {
        search_string.clear();
        search_string.append(exact_bcx);
        search_string.append(exact_bcy);
        search_string.append(exact_bcz);
        unsigned int endDistance;
        barcode.editd = static_cast<int>(
            edit_distance(Final_Barcode, search_string, endDistance, 8));
        barcode.umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index,
                              buindex.barcode_UMI_kmer_index.at(Final_Barcode));
        // FIXME(P1): barcodeX/Y/Z_unambiguous were computed and then ignored,
        // so this is unconditionally true whenever a whitelist hit exists.
        barcode.unambiguous = true;
        barcode.barcode = Final_Barcode;
    } else {
        bump(g_stats.whitelist_no_match);
    }

    return barcode;
}

// ------------------------------------------------------------------ visium --
inline Barcode get_visium_barcode(std::string& seq, const BarcodeSet* known_bcX,
                                  int flank_max_editd, const Pattern& search_pattern,
                                  const WhiteListMap* WHITELIST, BarcodeUMIindex& buindex) {
    Clock::time_point t1, t2;

    Barcode barcode;
    barcode.editd = -1000;
    barcode.flank_editd = 100;
    barcode.unambiguous = false;

    int n_eq = 0;
    const EdlibEqualityPair* eq = iupacEqualities(n_eq);
    EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, eq, n_eq};

    std::string search_string;
    search_string.reserve(60);
    std::vector<unsigned long> subpattern_lengths;
    subpattern_lengths.reserve(search_pattern.size());
    for (const auto& pair : search_pattern) {
        search_string += pair.second;
        subpattern_lengths.push_back(pair.second.length());
    }

    if (seq.length() < search_string.length()) {
        bump(g_stats.read_too_short);
        return barcode;
    }

    t1 = Clock::now();
    EdlibAlignResult result = edlibAlign(search_string.c_str(),
                                         static_cast<int>(search_string.length()),
                                         seq.c_str(), static_cast<int>(seq.length()),
                                         edlibConf);
    t2 = Clock::now();
    bump(g_stats.t_adapter_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (result.status != EDLIB_STATUS_OK || result.numLocations == 0) {
        bump(g_stats.adapter_not_found);
        edlibFreeAlignResult(result);
        return barcode;
    }

    barcode.flank_editd = result.editDistance;
    barcode.flank_start = result.startLocations[0];
    barcode.flank_end = result.endLocations[0];

    std::vector<int> read_to_subpatterns;
    const bool mapped = mapSubpatterns(result, barcode.flank_start,
                                       subpattern_lengths, read_to_subpatterns);
    edlibFreeAlignResult(result);
    if (!mapped) {
        bump(g_stats.adapter_not_found);
        return barcode;
    }

    const int bc_index = 1, umi_index = 2, edit_X = 5;

    std::string BarcodeX;
    std::vector<std::string> BarcodeVecStrX;

    t1 = Clock::now();
    std::string exact_bcx = seq.substr(
        static_cast<size_t>(read_to_subpatterns[bc_index]),
        static_cast<size_t>(read_to_subpatterns[bc_index + 1] - read_to_subpatterns[bc_index]));
    if (known_bcX->empty() || known_bcX->find(exact_bcx) != known_bcX->end()) {
        BarcodeX = exact_bcx;
        BarcodeVecStrX.push_back(BarcodeX);
    } else {
        std::unordered_map<std::string, int> candidate_votes;
        std::string bcx_ext;
        bcx_ext.reserve(exact_bcx.size() + 15);
        bcx_ext.append("CCGATCT");
        bcx_ext += exact_bcx;

        std::string kmer(visiumXKMER, '\0');
        for (size_t p = 0; p + visiumXKMER <= bcx_ext.size(); ++p) {
            std::memcpy(&kmer[0], bcx_ext.data() + p, visiumXKMER);
            auto it = buindex.barcodeX_kmer_index.find(kmer);
            if (it != buindex.barcodeX_kmer_index.end())
                for (const std::string& id : it->second) candidate_votes[id]++;
        }
        std::vector<std::pair<std::string, int>> cands(candidate_votes.begin(),
                                                       candidate_votes.end());
        std::sort(cands.begin(), cands.end(),
                  [](const std::pair<std::string, int>& a,
                     const std::pair<std::string, int>& b) { return a.second > b.second; });
        if (cands.size() > 10) cands.resize(10);
        for (const auto& c : cands) BarcodeVecStrX.push_back(c.first);
    }
    t2 = Clock::now();
    bump(g_stats.t_bc1_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    const int left_bound = std::max(read_to_subpatterns[bc_index] - 5, 0);
    const int max_length = static_cast<int>(search_pattern[bc_index].second.length()) + 2 * 5;
    const std::string barcode_seq = seq.substr(static_cast<size_t>(left_bound),
                                               static_cast<size_t>(max_length));

    if (BarcodeVecStrX.empty() && BarcodeX.empty()) {
        bump(g_stats.bc1_unresolved);
        unsigned int editDistance, endDistance;
        for (auto itr = known_bcX->begin(); itr != known_bcX->end(); ++itr) {
            search_string = *itr;
            editDistance = edit_distance(barcode_seq, search_string, endDistance, 3);
            if (static_cast<int>(editDistance) == barcode.editd) {
                barcode.unambiguous = false;
            } else if (static_cast<int>(editDistance) < barcode.editd && editDistance <= 3) {
                BarcodeVecStrX.clear();
                barcode.unambiguous = true;
                barcode.editd = static_cast<int>(editDistance);
                barcode.barcode = *itr;
                BarcodeVecStrX.push_back(barcode.barcode);
                if (editDistance == 0) break;
            }
        }
    }

    if (BarcodeVecStrX.size() == 1) exact_bcx = BarcodeVecStrX[0];

    t1 = Clock::now();
    std::string Final_Barcode;
    Final_Barcode.reserve(16);
    if (!WHITELIST->empty()) {
        wfa::WFAlignerGapAffine aligner(4, 5, 1, wfa::WFAligner::Alignment,
                                        wfa::WFAligner::MemoryHigh);
        for (const auto& fin : BarcodeVecStrX) {
            aligner.alignEndsFree(fin, 0, 0, barcode_seq, 1, 1);
            if (aligner.getAlignmentScore() > barcode.editd) {
                Final_Barcode = fin;
                barcode.editd = aligner.getAlignmentScore();
            }
        }
    }
    t2 = Clock::now();
    bump(g_stats.t_whitelist_us, std::chrono::duration_cast<microsec>(t2 - t1).count());

    if (!Final_Barcode.empty()) {
        unsigned int endDistance;
        barcode.editd = static_cast<int>(
            edit_distance(barcode_seq, Final_Barcode, endDistance, edit_X));
        barcode.umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index,
                              buindex.barcode_UMI_kmer_index.at(Final_Barcode));
        barcode.unambiguous = true;
        barcode.barcode = Final_Barcode;
    } else {
        bump(g_stats.whitelist_no_match);
    }

    return barcode;
}

// ----------------------------------------------------------- read dispatch --
inline std::vector<Barcode> big_barcode_search(const std::string& mode,
                                               std::string& sequence,
                                               const BarcodeSet& bcX, const BarcodeSet& bcY,
                                               const BarcodeSet& bcZ, int max_flank_editd,
                                               const Pattern& search_pattern,
                                               const WhiteListMap& whitelist,
                                               BarcodeUMIindex& buindex) {
    std::vector<Barcode> return_vec;
    Barcode result;

    if (mode == "magicseq") {
        result = get_magicseq_barcode(sequence, &bcX, &bcY, &bcZ, max_flank_editd,
                                      search_pattern, &whitelist, buindex);
        if (result.editd <= 2 && result.unambiguous) return_vec.push_back(result);
        if (!result.unambiguous) bump(g_stats.search_no_call);
        if (result.editd > 2) bump(g_stats.barcode_editd_high);
    } else if (mode == "visium") {
        result = get_visium_barcode(sequence, &bcX, max_flank_editd, search_pattern,
                                    &whitelist, buindex);
        if (result.editd <= 2 && result.unambiguous) return_vec.push_back(result);
        if (!result.unambiguous) bump(g_stats.search_no_call);
        if (result.editd > 2) bump(g_stats.barcode_editd_high);
    }

    // Mask what was found and look for further barcodes in the remainder.
    if (!return_vec.empty()) {
        std::string masked = sequence;
        bool masked_anything = false;
        for (const auto& b : return_vec) {
            const int flank_length = b.flank_end - b.flank_start;
            if (flank_length > 0 && b.flank_start >= 0 &&
                static_cast<size_t>(b.flank_start + flank_length) <= masked.size()) {
                masked.replace(static_cast<size_t>(b.flank_start),
                               static_cast<size_t>(flank_length),
                               std::string(static_cast<size_t>(flank_length), 'X'));
                masked_anything = true;
            }
        }
        // Guard, not a behaviour change: whenever the mask above is applied -
        // which is every well-formed hit - this recurses exactly as 1.x did.
        // What it no longer does is recurse when nothing could be masked; that
        // path searched an unchanged sequence, found the same hit again and
        // had no way to terminate.
        if (masked_anything) {
            std::vector<Barcode> more = big_barcode_search(mode, masked, bcX, bcY, bcZ,
                                                           max_flank_editd, search_pattern,
                                                           whitelist, buindex);
            return_vec.insert(return_vec.end(), more.begin(), more.end());
        }
    }
    return return_vec;
}

// One worker's slice of a batch.
inline void search_reads(const std::string& mode, std::vector<SearchResult>& reads,
                         const BarcodeSet& bcX, const BarcodeSet& bcY,
                         const BarcodeSet& bcZ, int flank_edit_distance,
                         const Pattern& search_pattern, const WhiteListMap& whitelist,
                         BarcodeUMIindex& buindex) {
    for (auto& r : reads) {
        r.vec_bc_for = big_barcode_search(mode, r.line, bcX, bcY, bcZ,
                                          flank_edit_distance, search_pattern,
                                          whitelist, buindex);
        r.rev_line = r.line;
        reverse_complement(r.rev_line);
        r.vec_bc_rev = big_barcode_search(mode, r.rev_line, bcX, bcY, bcZ,
                                          flank_edit_distance, search_pattern,
                                          whitelist, buindex);
        r.count = static_cast<int>(r.vec_bc_for.size() + r.vec_bc_rev.size());
        r.chimeric = !r.vec_bc_for.empty() && !r.vec_bc_rev.empty();
    }
}

}  // namespace prebrocoli

#endif  // PREBROCOLI_DEMUX_HPP
