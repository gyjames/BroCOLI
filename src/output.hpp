// Output layer.
//
// One abstraction covers both modes: results are written into *tables*, each
// table having a set of *columns*.
//   bulk : one table, one column per input BAM
//   sc   : one table per input BAM, one column per barcode
// Everything upstream only ever talks about (table, column) pairs, which is why
// the quantification code has a single implementation instead of the three
// near-identical copies of BroCOLI 1.x.
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.hpp"
#include "gtf.hpp"

namespace brocoli {

// ---------------------------------------------------------- column model ---
class ColumnSpace {
public:
    void init(Mode mode, const std::vector<std::string>& sample_names,
              const std::vector<std::set<std::string>>& barcodes_per_sample) {
        mode_ = mode;
        if (mode == Mode::Bulk) {
            tables_ = 1;
            table_tag_.push_back("");
            columns_.push_back(sample_names);
        } else {
            tables_ = static_cast<int>(sample_names.size());
            for (int i = 0; i < tables_; ++i) {
                table_tag_.push_back("_" + std::to_string(i));
                columns_.emplace_back(barcodes_per_sample[i].begin(),
                                      barcodes_per_sample[i].end());
                std::unordered_map<std::string, int> idx;
                idx.reserve(columns_[i].size() * 2);
                for (size_t j = 0; j < columns_[i].size(); ++j)
                    idx.emplace(columns_[i][j], static_cast<int>(j));
                bc_index_.push_back(std::move(idx));
            }
        }
    }

    int tables() const { return tables_; }
    int columns(int t) const { return static_cast<int>(columns_[t].size()); }
    const std::vector<std::string>& columnNames(int t) const { return columns_[t]; }
    const std::string& tag(int t) const { return table_tag_[t]; }

    int tableOf(int file) const { return mode_ == Mode::Bulk ? 0 : file; }

    int columnOf(int file, const std::string& bc) const {
        if (mode_ == Mode::Bulk) return file;
        auto it = bc_index_[file].find(bc);
        return it == bc_index_[file].end() ? -1 : it->second;
    }

private:
    Mode mode_ = Mode::Bulk;
    int tables_ = 1;
    std::vector<std::string> table_tag_;
    std::vector<std::vector<std::string>> columns_;
    std::vector<std::unordered_map<std::string, int>> bc_index_;
};

// ------------------------------------------------------------- row buffer --
// A sparse row: transcript (or gene) id plus (column, value) pairs. Keeping the
// intermediate sparse matters enormously in single-cell mode, where a dense row
// would be tens of thousands of mostly-zero cells.
struct SparseRow {
    std::string name;
    std::string gene;              // empty for gene rows
    std::vector<std::pair<int, double>> vals;
};

class TableSink {
public:
    TableSink(const std::string& path) : path_(path) {
        out_.open(path, std::ios::trunc);
        if (!out_) throw std::runtime_error("cannot open " + path);
    }
    const std::string& path() const { return path_; }

    void write(const std::vector<SparseRow>& rows) {
        if (rows.empty()) return;
        std::string buf;
        buf.reserve(rows.size() * 64);
        for (const auto& r : rows) {
            if (r.vals.empty()) continue;
            buf += r.name;
            buf += '\t';
            buf += r.gene.empty() ? "." : r.gene;
            for (const auto& v : r.vals) {
                buf += '\t';
                buf += std::to_string(v.first);
                buf += ':';
                buf += std::to_string(v.second);
            }
            buf += '\n';
        }
        if (buf.empty()) return;
        std::lock_guard<std::mutex> lk(mu_);
        out_ << buf;
    }

    void close() { std::lock_guard<std::mutex> lk(mu_); out_.close(); }

private:
    std::string path_;
    std::ofstream out_;
    std::mutex mu_;
};

class TextSink {
public:
    TextSink(const std::string& path, const std::string& header) : path_(path) {
        out_.open(path, std::ios::trunc);
        if (!out_) throw std::runtime_error("cannot open " + path);
        if (!header.empty()) out_ << header << '\n';
    }
    void write(const std::string& buf) {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> lk(mu_);
        out_ << buf;
    }
    bool isOpen() const { return out_.is_open(); }
    void close() { std::lock_guard<std::mutex> lk(mu_); out_.close(); }

private:
    std::string path_;
    std::ofstream out_;
    std::mutex mu_;
};

// ------------------------------------------------------------ gtf records --
inline void appendAnnotatedGtf(std::string& buf, const std::string& chrom,
                               const std::string& txKey, const IntervalVec& exons,
                               const std::string& strand) {
    if (exons.empty()) return;
    std::string gene, tx;
    splitTxKey(txKey, gene, tx);
    const std::string st = strand.empty() ? "." : strand;

    buf += chrom + "\tannotated_isoform\ttranscript\t" +
           std::to_string(exons.front()[0]) + '\t' + std::to_string(exons.back()[1]) +
           "\t.\t" + st + "\t.\tgene_id \"" + gene + "\"; transcript_id \"" + tx + "\";\n";
    for (size_t i = 0; i < exons.size(); ++i)
        buf += chrom + "\tannotated_isoform\texon\t" +
               std::to_string(exons[i][0]) + '\t' + std::to_string(exons[i][1]) +
               "\t.\t" + st + "\t.\tgene_id \"" + gene + "\"; transcript_id \"" + tx +
               "\"; exon_number \"" + std::to_string(i) + "\";\n";
}

inline void appendNovelGtf(std::string& buf, const std::string& chrom,
                           const std::string& txName, const std::string& gene,
                           const Interval& span, const IntervalVec& sjs,
                           const std::string& strand) {
    if (sjs.empty()) return;
    const std::string st = strand.empty() ? "." : strand;
    buf += chrom + "\tnovel_isoform\ttranscript\t" + std::to_string(span[0]) + '\t' +
           std::to_string(span[1]) + "\t.\t" + st + "\t.\tgene_id \"" + gene +
           "\"; transcript_id \"" + txName + "\";\n";

    for (size_t i = 0; i <= sjs.size(); ++i) {
        const int b = (i == 0) ? span[0] : sjs[i - 1][1] + 1;
        const int e = (i == sjs.size()) ? span[1] : sjs[i][0] - 1;
        buf += chrom + "\tnovel_isoform\texon\t" + std::to_string(b) + '\t' +
               std::to_string(e) + "\t.\t" + st + "\t.\tgene_id \"" + gene +
               "\"; transcript_id \"" + txName + "\"; exon_number \"" +
               std::to_string(i) + "\";\n";
    }
}

// --------------------------------------------------------- final rewrite ---
// The per-group writers can emit the same transcript from more than one group
// (a gene straddling a cluster boundary, or a single-exon model that also
// collected spliced reads), so the sparse intermediate is folded once at the
// end. This also converts to the requested output format.
struct RewriteOptions {
    bool mtx = false;
    int  min_count = 0;
    bool with_gene_column = true;
};

inline void rewriteTable(const std::string& sparse_path, const std::string& out_prefix,
                         const std::vector<std::string>& column_names,
                         const std::string& id_header, const RewriteOptions& ro) {
    std::unordered_map<std::string, size_t> index;
    std::vector<std::string> names, genes;
    std::vector<std::map<int, double>> data;

    std::ifstream in(sparse_path);
    if (!in) {
        std::cerr << "[BroCOLI] cannot re-read " << sparse_path << "\n";
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        size_t p = line.find('\t');
        if (p == std::string::npos) continue;
        std::string name = line.substr(0, p);
        size_t q = line.find('\t', p + 1);
        std::string gene = (q == std::string::npos) ? "." : line.substr(p + 1, q - p - 1);
        if (q == std::string::npos) q = line.size();

        auto it = index.find(name);
        size_t row;
        if (it == index.end()) {
            row = names.size();
            index.emplace(name, row);
            names.push_back(std::move(name));
            genes.push_back(std::move(gene));
            data.emplace_back();
        } else {
            row = it->second;
        }

        size_t pos = q;
        while (pos < line.size()) {
            const size_t start = pos + 1;
            if (start >= line.size()) break;
            size_t tab = line.find('\t', start);
            if (tab == std::string::npos) tab = line.size();
            const size_t colon = line.find(':', start);
            if (colon != std::string::npos && colon < tab) {
                const int c = std::atoi(line.c_str() + start);
                const double v = std::atof(line.c_str() + colon + 1);
                data[row][c] += v;
            }
            pos = tab;
        }
    }
    in.close();

    auto keep = [&](double v) { return v >= ro.min_count && v > 0.0; };

    if (!ro.mtx) {
        std::ofstream out(out_prefix + ".txt", std::ios::trunc);
        std::string buf;
        buf.reserve(1 << 20);
        buf += id_header;
        if (ro.with_gene_column) buf += "\tgene_id";
        for (const auto& c : column_names) { buf += '\t'; buf += c; }
        buf += '\n';
        out << buf;

        const size_t ncol = column_names.size();
        for (size_t r = 0; r < names.size(); ++r) {
            buf.clear();
            buf += names[r];
            if (ro.with_gene_column) { buf += '\t'; buf += genes[r]; }
            for (size_t c = 0; c < ncol; ++c) {
                auto it = data[r].find(static_cast<int>(c));
                double v = (it == data[r].end()) ? 0.0 : it->second;
                if (!keep(v)) v = 0.0;
                buf += '\t';
                buf += std::to_string(static_cast<long long>(std::llround(v)));
            }
            buf += '\n';
            out << buf;
        }
        out.close();
        return;
    }

    // MatrixMarket triplet output (features x columns), 1-based indices.
    {
        std::ofstream f(out_prefix + "_features.tsv", std::ios::trunc);
        for (size_t r = 0; r < names.size(); ++r) {
            f << names[r];
            if (ro.with_gene_column) f << '\t' << genes[r];
            f << '\n';
        }
    }
    {
        std::ofstream b(out_prefix + "_barcodes.tsv", std::ios::trunc);
        for (const auto& c : column_names) b << c << '\n';
    }
    long long nnz = 0;
    for (const auto& row : data)
        for (const auto& kv : row) if (keep(kv.second)) ++nnz;

    std::ofstream m(out_prefix + "_matrix.mtx", std::ios::trunc);
    m << "%%MatrixMarket matrix coordinate integer general\n";
    m << names.size() << ' ' << column_names.size() << ' ' << nnz << '\n';
    std::string buf;
    buf.reserve(1 << 20);
    for (size_t r = 0; r < data.size(); ++r) {
        for (const auto& kv : data[r]) {
            if (!keep(kv.second)) continue;
            buf += std::to_string(r + 1);
            buf += ' ';
            buf += std::to_string(kv.first + 1);
            buf += ' ';
            buf += std::to_string(static_cast<long long>(std::llround(kv.second)));
            buf += '\n';
            if (buf.size() > (1u << 20)) { m << buf; buf.clear(); }
        }
    }
    m << buf;
}

}  // namespace brocoli
