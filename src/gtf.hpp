// GTF parsing.
//
// Difference to BroCOLI 1.x: the old parser was a stateful streaming scanner
// that assumed every transcript's exons were contiguous in the file and that a
// descending exon order meant a minus-strand transcript. Here we bucket exons
// per (chromosome, gene|transcript) and sort at the end, which is order
// independent and therefore safe with any GENCODE/Ensembl/RefSeq flavour. The
// GTF is also only materialised once (1.x built a std::map copy and then an
// unordered_map copy of the same thing).
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.hpp"

namespace brocoli {

struct ChromAnno {
    // key = "gene_id|transcript_id"
    std::unordered_map<std::string, IntervalVec> tx_exons;
    std::unordered_map<std::string, IntervalVec> tx_sjs;    // multi-exon only
    std::unordered_map<std::string, Interval>    tx_span;
    std::unordered_map<std::string, std::string> tx_strand;

    std::unordered_map<std::string, Interval>    gene_span;
    std::unordered_map<std::string, std::string> gene_strand;
    std::unordered_map<std::string, std::vector<std::string>> gene2tx;
};

struct Annotation {
    std::unordered_map<std::string, ChromAnno> chrom;

    const ChromAnno* find(const std::string& c) const {
        auto it = chrom.find(c);
        return it == chrom.end() ? nullptr : &it->second;
    }
    size_t transcriptCount() const {
        size_t n = 0;
        for (const auto& kv : chrom) n += kv.second.tx_exons.size();
        return n;
    }
};

inline std::string extractGtfAttr(const std::string& attrs, const std::string& key) {
    const std::string needle = key + " \"";
    size_t pos = attrs.find(needle);
    if (pos == std::string::npos) return {};
    pos += needle.size();
    const size_t end = attrs.find('"', pos);
    return end == std::string::npos ? std::string{} : attrs.substr(pos, end - pos);
}

inline Annotation loadGtf(const std::string& path) {
    Annotation anno;
    if (path.empty()) {
        std::cerr << "[BroCOLI] no GTF supplied - running in annotation-free mode\n";
        return anno;
    }

    std::ifstream in(path);
    if (!in) {
        std::cerr << "[BroCOLI] cannot open GTF: " << path << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cout << "[BroCOLI] reading GTF " << path << "\n";

    std::string line, chr, feature, attrs, strand;
    long long nexon = 0;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Split on tabs without building a stringstream: the attribute column
        // is everything after the 8th tab.
        const char* p = line.c_str();
        const char* col[9];
        size_t ncol = 0;
        col[ncol++] = p;
        for (const char* q = p; *q && ncol < 9; ++q)
            if (*q == '\t') col[ncol++] = q + 1;
        if (ncol < 9) continue;

        auto field = [&](size_t i) {
            const char* s = col[i];
            const char* e = (i + 1 < ncol) ? col[i + 1] - 1 : line.c_str() + line.size();
            return std::string(s, e - s);
        };

        feature = field(2);
        if (feature != "exon") continue;

        chr    = field(0);
        strand = field(6);
        attrs.assign(col[8], line.c_str() + line.size() - col[8]);

        const std::string gene = extractGtfAttr(attrs, "gene_id");
        const std::string tx   = extractGtfAttr(attrs, "transcript_id");
        if (gene.empty() || tx.empty()) continue;

        const std::string key = gene + "|" + tx;
        ChromAnno& ca = anno.chrom[chr];
        ca.tx_exons[key].push_back(Interval{std::atoi(field(3).c_str()),
                                            std::atoi(field(4).c_str())});
        ca.tx_strand[key] = strand;
        ++nexon;
    }

    // ------------------------------------------------ derive SJs and genes --
    for (auto& ckv : anno.chrom) {
        ChromAnno& ca = ckv.second;

        for (auto& tkv : ca.tx_exons) {
            IntervalVec& ex = tkv.second;
            std::sort(ex.begin(), ex.end(),
                      [](const Interval& a, const Interval& b) { return a[0] < b[0]; });

            ca.tx_span[tkv.first] = Interval{ex.front()[0], ex.back()[1]};

            if (ex.size() > 1) {
                IntervalVec sjs;
                sjs.reserve(ex.size() - 1);
                for (size_t i = 0; i + 1 < ex.size(); ++i)
                    sjs.push_back(Interval{ex[i][1] + 1, ex[i + 1][0] - 1});
                ca.tx_sjs[tkv.first] = std::move(sjs);
            }

            std::string gene, txid;
            splitTxKey(tkv.first, gene, txid);
            ca.gene2tx[gene].push_back(tkv.first);

            auto git = ca.gene_span.find(gene);
            const Interval& span = ca.tx_span[tkv.first];
            if (git == ca.gene_span.end()) {
                ca.gene_span[gene]   = span;
                ca.gene_strand[gene] = ca.tx_strand[tkv.first];
            } else {
                git->second[0] = std::min(git->second[0], span[0]);
                git->second[1] = std::max(git->second[1], span[1]);
            }
        }
    }

    size_t ngene = 0;
    for (const auto& kv : anno.chrom) ngene += kv.second.gene_span.size();
    std::cout << "[BroCOLI] GTF: " << nexon << " exons, "
              << anno.transcriptCount() << " transcripts, " << ngene << " genes\n";
    return anno;
}

// ------------------------------------------------------------- per group ---
// The subset of the annotation that overlaps one read cluster.
struct GroupAnno {
    std::unordered_map<std::string, Interval>    genes;
    std::unordered_map<std::string, IntervalVec> me_exons;  // multi-exon transcripts
    std::unordered_map<std::string, IntervalVec> me_sjs;
    std::unordered_map<std::string, Interval>    me_span;
    std::unordered_map<std::string, Interval>    se_tx;     // single-exon transcripts
    // gene -> merged exon blocks, used for single-exon read assignment
    std::unordered_map<std::string, IntervalVec> gene_exons;
};

inline GroupAnno sliceAnnotation(const ChromAnno& ca, const Interval& span) {
    GroupAnno g;
    for (const auto& kv : ca.gene_span) {
        if (kv.second[1] < span[0] || span[1] < kv.second[0]) continue;
        g.genes.emplace(kv.first, kv.second);
    }
    if (g.genes.empty()) return g;

    g.me_exons.reserve(g.genes.size() * 4);
    g.me_sjs.reserve(g.genes.size() * 4);
    g.me_span.reserve(g.genes.size() * 4);
    g.se_tx.reserve(g.genes.size() * 2);
    g.gene_exons.reserve(g.genes.size());

    for (const auto& gk : g.genes) {
        auto it = ca.gene2tx.find(gk.first);
        if (it == ca.gene2tx.end()) continue;

        IntervalVec allExons;
        for (const std::string& tx : it->second) {
            auto ex = ca.tx_exons.find(tx);
            if (ex == ca.tx_exons.end() || ex->second.empty()) continue;
            allExons.insert(allExons.end(), ex->second.begin(), ex->second.end());

            if (ex->second.size() > 1) {
                g.me_exons.emplace(tx, ex->second);
                auto sj = ca.tx_sjs.find(tx);
                if (sj != ca.tx_sjs.end()) g.me_sjs.emplace(tx, sj->second);
                auto sp = ca.tx_span.find(tx);
                if (sp != ca.tx_span.end()) g.me_span.emplace(tx, sp->second);
            } else {
                g.se_tx.emplace(tx, ex->second[0]);
            }
        }
        g.gene_exons.emplace(gk.first, mergeIntervals(std::move(allExons)));
    }
    return g;
}

}  // namespace brocoli
