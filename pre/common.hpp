// preBroCOLI 2.0 - options, statistics and small shared types.
#ifndef PREBROCOLI_COMMON_HPP
#define PREBROCOLI_COMMON_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace prebrocoli {

inline const char* kVersion = "2.0.0";

// --------------------------------------------------------------- options ---
struct Options {
    std::string mode;              // magicseq | 10x3v3 | visium
    std::string reads;             // "-" for stdin; .fastq or .fastq.gz
    std::string whitelist;
    std::string barcodeX, barcodeY, barcodeZ;

    std::string outdir = ".";
    std::string prefix = "preBroCOLI";

    int  threads          = 1;
    int  batch            = 2000;   // reads per thread per batch
    int  flank_editd      = 20;     // overridden per chemistry unless given
    bool flank_editd_set  = false;

    bool trim_barcodes    = true;   // -i
    bool split_by_barcode = false;  // -s
    bool mark_chimeric    = false;  // -c
    bool keep_unmatched   = false;  // --keep_unmatched
};

// ------------------------------------------------------------- statistics --
//
// Every counter here used to be a bare int incremented from worker threads,
// i.e. a data race: the numbers printed were too low and differed between
// runs. They are atomics now.
//
// Two levels are counted, and they are NOT the same thing:
//
//   * read level    - one count per input read (reads_in, reads_assigned,
//                     reads_unmatched, reads_chimeric, reads_written).
//                     reads_in == reads_assigned + reads_unmatched always;
//                     this invariant is asserted in the tests.
//   * attempt level - one count per barcode search. Each read is searched on
//                     both strands and searches recurse after a hit, so these
//                     do not sum to reads_in. They tell you *where* searches
//                     die, not how many reads died.
//
// The t_* timers are summed across threads, so they are CPU time, not wall
// time: 32 threads for 10 minutes prints as ~19000 s. Labelled as such.
struct Stats {
    // read level
    std::atomic<long long> reads_in{0};
    std::atomic<long long> reads_assigned{0};
    std::atomic<long long> reads_unmatched{0};
    std::atomic<long long> reads_chimeric{0};
    std::atomic<long long> reads_multi_barcode{0};
    std::atomic<long long> reads_written{0};

    // attempt level - why a barcode search failed
    std::atomic<long long> read_too_short{0};      // was noLen
    std::atomic<long long> adapter_not_found{0};   // was noAlign
    std::atomic<long long> bc1_unresolved{0};      // was noX
    std::atomic<long long> bc2_unresolved{0};      // was noY
    std::atomic<long long> bc3_unresolved{0};      // was noZ
    std::atomic<long long> whitelist_no_match{0};  // was noFinal
    std::atomic<long long> search_no_call{0};      // was Ambiguous - counts every search
                                                   // that produced no confident barcode,
                                                   // not only genuine ties
    std::atomic<long long> barcode_editd_high{0};  // was total6 (editd > 2)

    // UMI resolution
    std::atomic<long long> umi_exact{0};           // was erfen
    std::atomic<long long> umi_corrected{0};       // was bitwei
    std::atomic<long long> umi_uncorrected{0};     // was chushi

    // CPU time, microseconds, summed over threads
    std::atomic<long long> t_adapter_us{0};        // was total_First
    std::atomic<long long> t_bc1_us{0};            // was total_BarcodeX
    std::atomic<long long> t_bc2_us{0};            // was total_BarcodeY
    std::atomic<long long> t_bc3_us{0};            // was total_BarcodeZ
    std::atomic<long long> t_whitelist_us{0};      // was total_Final
};

inline Stats g_stats;

inline void bump(std::atomic<long long>& c, long long by = 1) {
    c.fetch_add(by, std::memory_order_relaxed);
}

// ------------------------------------------------------------ small types --
struct Barcode {
    std::string barcode;
    std::string umi;
    int editd = 0;
    int flank_editd = 0;
    int flank_start = 0;
    int flank_end = 0;
    bool unambiguous = false;
};

struct SearchResult {
    std::string read_id;
    std::string qual_scores;
    std::string line;
    std::string rev_line;
    std::vector<Barcode> vec_bc_for;
    std::vector<Barcode> vec_bc_rev;
    int count = 0;
    bool chimeric = false;
};

// ---------------------------------------------------------------- utility --
inline bool fileExists(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    if (!S_ISREG(st.st_mode)) {
        std::cerr << "[preBroCOLI] not a regular file: " << path << "\n";
        return false;
    }
    if (access(path.c_str(), R_OK) != 0) {
        std::cerr << "[preBroCOLI] file exists but is not readable: " << path << "\n";
        return false;
    }
    return true;
}

inline bool makeDir(const std::string& path) {
    if (path.empty() || path == ".") return true;
    struct stat st;
    if (stat(path.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    if (mkdir(path.c_str(), 0755) == 0) return true;
    std::cerr << "[preBroCOLI] cannot create " << path << ": "
              << std::strerror(errno) << "\n";
    return false;
}

inline std::string joinPath(const std::string& a, const std::string& b) {
    if (a.empty() || a == ".") return b;
    return a.back() == '/' ? a + b : a + "/" + b;
}

static const std::array<unsigned char, 256> lut = [] {
    std::array<unsigned char, 256> t;
    t.fill(0xFF);
    t[static_cast<unsigned char>('A')] = 0; t[static_cast<unsigned char>('a')] = 0;
    t[static_cast<unsigned char>('C')] = 1; t[static_cast<unsigned char>('c')] = 1;
    t[static_cast<unsigned char>('G')] = 2; t[static_cast<unsigned char>('g')] = 2;
    t[static_cast<unsigned char>('T')] = 3; t[static_cast<unsigned char>('t')] = 3;
    return t;
}();

inline bool encode_dna_2bit_n(const char* seq, int len, uint32_t& out) {
    out = 0;
    for (int i = 0; i < len; ++i) {
        const unsigned char v = lut[static_cast<unsigned char>(seq[i])];
        if (v == 0xFF) return false;
        out = (out << 2) | v;
    }
    return true;
}

inline char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default:  return 'N';
    }
}

inline void reverse_complement(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
    std::transform(seq.begin(), seq.end(), seq.begin(),
                   [](char c) { return complement(c); });
}

}  // namespace prebrocoli

#endif  // PREBROCOLI_COMMON_HPP
