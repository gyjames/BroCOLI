// preBroCOLI 2.0 - whitelist / barcode loading and k-mer index construction.
//
// Behaviour is carried over unchanged from 1.x so that P0 output matches the
// old binary byte for byte. The differences are in error handling only:
//   * the magicseq branch used to test file_exists(barcodePathX) three times,
//     so a missing Y or Z file was reported as "0 known barcodesY" instead of
//     "file not found";
//   * unimplemented chemistries used to print a note and carry on with empty
//     tables, producing zero barcodes and exit status 0. They exit(1) now.
#ifndef PREBROCOLI_WHITELIST_HPP
#define PREBROCOLI_WHITELIST_HPP

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "common.hpp"

namespace prebrocoli {

constexpr int magicseqXKMER  = 3;
constexpr int magicseqALLKMER = 9;
constexpr int visiumXKMER    = 11;
constexpr int visiumALLKMER  = 11;

struct UmiIndex {
    std::vector<uint32_t> valid_umis;
    std::vector<std::vector<int>> kmer_lut;
    int k = 5;
    int umi_len = 12;

    static std::string decode32(uint32_t val, int len) {
        std::string res(static_cast<size_t>(len), 'N');
        static const char dec[] = "ACGT";
        for (int i = len - 1; i >= 0; --i) {
            res[static_cast<size_t>(i)] = dec[val & 3];
            val >>= 2;
        }
        return res;
    }

    void build(const std::vector<uint32_t>& umis, int k_val = 5, int len_val = 12) {
        k = k_val;
        umi_len = len_val;
        valid_umis = umis;

        const int num_kmers = 1 << (2 * k);
        kmer_lut.clear();
        kmer_lut.resize(static_cast<size_t>(num_kmers));

        const uint32_t mask = (1U << (2 * k)) - 1;
        for (size_t i = 0; i < valid_umis.size(); ++i) {
            const uint32_t u = valid_umis[i];
            for (int p = 0; p <= umi_len - k; ++p) {
                const int shift = 2 * (umi_len - p - k);
                const uint32_t kmer_val = (u >> shift) & mask;
                kmer_lut[kmer_val].push_back(static_cast<int>(i));
            }
        }
    }
};

struct BarcodeUMIindex {
    std::unordered_map<std::string, std::vector<std::string>> barcodeX_kmer_index;
    std::unordered_map<std::string, std::vector<std::string>> barcodeY_kmer_index;
    std::unordered_map<std::string, std::vector<std::string>> barcodeZ_kmer_index;
    std::unordered_map<std::string, std::vector<std::string>> whiteList_kmer_index;
    std::unordered_map<std::string, UmiIndex> barcode_UMI_kmer_index;
};

using WhiteListMap = std::unordered_map<std::string, std::vector<uint32_t>>;
using BarcodeSet   = std::unordered_set<std::string>;

// ---------------------------------------------------------------- loading --
namespace detail {

inline void dedupWhitelistUmis(WhiteListMap& wl, const char* tag) {
    size_t total = 0;
    for (auto& pair : wl) {
        std::vector<uint32_t>& v = pair.second;
        if (!v.empty()) {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            v.shrink_to_fit();
        }
        total += v.size();
    }
    std::cout << "[preBroCOLI] " << tag << " whitelist barcodes: " << wl.size() << "\n";
    std::cout << "[preBroCOLI] unique UMIs loaded: " << total << "\n";
    if (wl.empty()) {
        std::cerr << "[preBroCOLI] the whitelist file provides no barcodes\n";
        std::exit(1);
    }
}

// One column of barcodes, one per line (extra columns ignored).
inline void loadBarcodeColumn(const std::string& path, size_t min_len,
                              bool strip_suffix, BarcodeSet& out,
                              const char* label) {
    if (!fileExists(path)) {
        std::cerr << "[preBroCOLI] " << label << " file not found: " << path << "\n";
        std::exit(1);
    }
    std::ifstream file(path);
    std::cout << "[preBroCOLI] " << label << " from: " << path << "\n";
    std::string line, token;
    while (std::getline(file, line)) {
        std::istringstream ls(line);
        while (ls >> token) {
            if (token.size() < min_len) continue;
            if (strip_suffix) {
                const size_t pos = token.rfind('-');
                if (pos != std::string::npos) {
                    bool all_digit = true;
                    for (size_t i = pos + 1; i < token.size(); ++i)
                        if (!isdigit(static_cast<unsigned char>(token[i]))) { all_digit = false; break; }
                    if (all_digit) token = token.substr(0, pos);
                }
                if (token.size() < min_len) continue;
            }
            out.insert(token);
        }
    }
    std::cout << "[preBroCOLI] known " << label << ": " << out.size() << "\n";
    if (out.empty()) {
        std::cerr << "[preBroCOLI] " << label << " list is empty\n";
        std::exit(1);
    }
}

}  // namespace detail

inline void loadMagicseq(const Options& opt, WhiteListMap& wl,
                         BarcodeSet& bcX, BarcodeSet& bcY, BarcodeSet& bcZ) {
    if (fileExists(opt.whitelist)) {
        std::ifstream file(opt.whitelist);
        std::cout << "[preBroCOLI] [magicseq] barcodes X/Y/Z and UMIs from: "
                  << opt.whitelist << "\n";
        std::string line;
        uint32_t umi_val = 0;
        while (std::getline(file, line)) {
            if (line.size() < 24) continue;
            bcX.emplace(line, 0, 8);
            bcY.emplace(line, 8, 8);
            bcZ.emplace(line, 16, 8);
            std::vector<uint32_t>& umi_list = wl[line.substr(0, 24)];
            if (line.size() > 36 && encode_dna_2bit_n(line.data() + 25, 12, umi_val))
                umi_list.push_back(umi_val);
        }
        detail::dedupWhitelistUmis(wl, "[magicseq]");
        return;
    }

    std::cout << "[preBroCOLI] no whitelist file (" << opt.whitelist
              << "), falling back to per-position barcode lists\n";
    if (opt.barcodeX.empty() || opt.barcodeY.empty() || opt.barcodeZ.empty()) {
        std::cerr << "[preBroCOLI] magicseq needs either -w <whitelist> or all of "
                     "-x/-y/-z\n";
        std::exit(1);
    }
    detail::loadBarcodeColumn(opt.barcodeX, 8, false, bcX, "barcodeX");
    detail::loadBarcodeColumn(opt.barcodeY, 8, false, bcY, "barcodeY");
    detail::loadBarcodeColumn(opt.barcodeZ, 8, false, bcZ, "barcodeZ");
}

inline void loadSingleBarcode(const Options& opt, WhiteListMap& wl,
                              BarcodeSet& bcX, const char* tag) {
    if (fileExists(opt.whitelist)) {
        std::ifstream file(opt.whitelist);
        std::cout << "[preBroCOLI] " << tag << " barcodes and UMIs from: "
                  << opt.whitelist << "\n";
        std::string line;
        uint32_t umi_val = 0;
        while (std::getline(file, line)) {
            if (line.size() < 16) continue;
            const std::string bc = line.substr(0, 16);
            bcX.emplace(bc);
            std::vector<uint32_t>& umi_list = wl[bc];
            if (line.size() > 28 && encode_dna_2bit_n(line.data() + 17, 12, umi_val))
                umi_list.push_back(umi_val);
        }
        detail::dedupWhitelistUmis(wl, tag);
        return;
    }

    std::cout << "[preBroCOLI] no whitelist file (" << opt.whitelist
              << "), falling back to a barcode list\n";
    if (opt.barcodeX.empty()) {
        std::cerr << "[preBroCOLI] " << tag << " needs either -w <whitelist> or -x <barcodes>\n";
        std::exit(1);
    }
    detail::loadBarcodeColumn(opt.barcodeX, 16, true, bcX, "barcode");
}

// ---------------------------------------------------------------- indexes --
inline void buildMagicseqIndex(const WhiteListMap& wl, const BarcodeSet& bcX,
                               const BarcodeSet& bcY, const BarcodeSet& bcZ,
                               BarcodeUMIindex& ix) {
    std::string kmer;
    kmer.resize(magicseqXKMER);
    std::string bc_ext;
    bc_ext.reserve(16);

    std::cout << "[preBroCOLI] building [magicseq] k-mer indexes for X/Y/Z\n";

    for (const auto& bc : bcX) {
        bc_ext.clear();
        bc_ext.append("T"); bc_ext += bc; bc_ext.append("C");
        for (size_t p = 0; p + magicseqXKMER <= bc_ext.size(); ++p) {
            std::memcpy(&kmer[0], bc_ext.data() + p, magicseqXKMER);
            ix.barcodeX_kmer_index[kmer].emplace_back(bc);
        }
    }
    for (const auto& bc : bcY) {
        bc_ext.clear();
        bc_ext.append("A"); bc_ext += bc; bc_ext.append("T");
        for (size_t p = 0; p + magicseqXKMER <= bc_ext.size(); ++p) {
            std::memcpy(&kmer[0], bc_ext.data() + p, magicseqXKMER);
            ix.barcodeY_kmer_index[kmer].emplace_back(bc);
        }
    }
    for (const auto& bc : bcZ) {
        for (size_t p = 0; p + magicseqXKMER <= bc.size(); ++p) {
            std::memcpy(&kmer[0], bc.data() + p, magicseqXKMER);
            ix.barcodeZ_kmer_index[kmer].emplace_back(bc);
        }
    }

    kmer.resize(magicseqALLKMER);
    for (const auto& w : wl) {
        for (size_t p = 0; p + magicseqALLKMER <= w.first.size(); ++p) {
            std::memcpy(&kmer[0], w.first.data() + p, magicseqALLKMER);
            ix.whiteList_kmer_index[kmer].emplace_back(w.first);
        }
    }

    std::cout << "[preBroCOLI] building UMI k-mer indexes for " << wl.size()
              << " barcodes\n";
    for (const auto& pair : wl)
        ix.barcode_UMI_kmer_index[pair.first].build(pair.second, 4, 12);
    std::cout << "[preBroCOLI] index build complete\n";
}

// 10x3v3 and visium share the layout; only the primer constant differs, and
// that lives in the search pattern, not here.
inline void buildSingleBarcodeIndex(const WhiteListMap& wl, const BarcodeSet& bcX,
                                    const std::string& primer_tail,
                                    BarcodeUMIindex& ix, const char* tag) {
    std::string kmer;
    kmer.resize(visiumXKMER);
    std::string bc_ext;
    bc_ext.reserve(32);

    std::cout << "[preBroCOLI] building " << tag << " k-mer index for barcodes\n";
    for (const auto& bc : bcX) {
        bc_ext.clear();
        bc_ext.append(primer_tail);
        bc_ext += bc;
        for (size_t p = 0; p + visiumXKMER <= bc_ext.size(); ++p) {
            std::memcpy(&kmer[0], bc_ext.data() + p, visiumXKMER);
            ix.barcodeX_kmer_index[kmer].emplace_back(bc);
        }
    }

    kmer.resize(visiumALLKMER);
    for (const auto& w : wl) {
        for (size_t p = 0; p + visiumALLKMER <= w.first.size(); ++p) {
            std::memcpy(&kmer[0], w.first.data() + p, visiumALLKMER);
            ix.whiteList_kmer_index[kmer].emplace_back(w.first);
        }
    }

    std::cout << "[preBroCOLI] building UMI k-mer indexes for " << wl.size()
              << " barcodes\n";
    for (const auto& pair : wl)
        ix.barcode_UMI_kmer_index[pair.first].build(pair.second, 4, 12);
    std::cout << "[preBroCOLI] index build complete\n";
}

// --------------------------------------------------------- search pattern --
using Pattern = std::vector<std::pair<std::string, std::string>>;

inline Pattern searchPattern(const std::string& mode) {
    if (mode == "magicseq") {
        return {{"primer",   "CTACACGACGCTCTTCCGATCT"},
                {"BCX",      std::string(8, '?')},
                {"LinkerXY", "CAGTCATGTCATGAGCTA"},
                {"BCY",      std::string(8, '?')},
                {"LinkerYZ", "TGATGCGACACTGATCGA"},
                {"BCZ",      std::string(8, '?')},
                {"UMI",      std::string(12, '?')},
                {"polyT",    "TTTTTTTTTT"}};
    }
    // 10x3v3 and visium
    return {{"primer", "CTACACGACGCTCTTCCGATCT"},
            {"BC",     std::string(16, '?')},
            {"UMI",    std::string(12, '?')},
            {"polyT",  "TTTTTTTTTT"}};
}

inline void loadChemistry(const Options& opt, WhiteListMap& wl, BarcodeSet& bcX,
                          BarcodeSet& bcY, BarcodeSet& bcZ, BarcodeUMIindex& ix) {
    if (opt.mode == "magicseq") {
        std::cout << "[preBroCOLI] chemistry: MAGIC-seq\n";
        loadMagicseq(opt, wl, bcX, bcY, bcZ);
        ix.whiteList_kmer_index.reserve(wl.size());
        buildMagicseqIndex(wl, bcX, bcY, bcZ, ix);
    } else if (opt.mode == "visium") {
        std::cout << "[preBroCOLI] chemistry: Visium\n";
        loadSingleBarcode(opt, wl, bcX, "[visium]");
        ix.whiteList_kmer_index.reserve(wl.size());
        buildSingleBarcodeIndex(wl, bcX, "CCGATCT", ix, "[visium]");
    } else if (opt.mode == "10x3v3") {
        // The demux path for this chemistry is still empty (P2 completes it).
        // 1.x loaded the whitelist here and then silently produced no barcodes
        // at all; failing loudly is better than that.
        std::cerr << "[preBroCOLI] chemistry 10x3v3 is not implemented yet "
                     "(the demux branch is empty); use --mode visium if the "
                     "layout matches, or wait for the next release\n";
        std::exit(1);
    } else {
        std::cerr << "[preBroCOLI] unknown chemistry '" << opt.mode
                  << "'. Supported: magicseq, visium (10x3v3, HD, stereoseq "
                     "not implemented yet)\n";
        std::exit(1);
    }
}

}  // namespace prebrocoli

#endif  // PREBROCOLI_WHITELIST_HPP
