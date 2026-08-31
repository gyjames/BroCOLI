// Stage 1: coordinate-sorted BAM -> intermediate per-read records.
//
// Parallelism comes from the BAM index: the genome is cut into windows and one
// worker owns each window. Compared with the byte-offset chunking of BroCOLI
// 1.x this removes findNextLineStart(), the "stuck position" watchdog and the
// possibility of splitting a record in half. A read belongs to the window that
// contains its start position, so no read is emitted twice even though the
// iterator returns everything that *overlaps* the window.
//
// Stage 1 no longer forms read clusters: clustering happens once, globally,
// while merging (see group.hpp). That removes the window-boundary special case
// entirely.
#pragma once

#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "common.hpp"
#include "record.hpp"
#include "tpool.hpp"

namespace brocoli {

// ------------------------------------------------------------------ fasta --
using FastaRef = std::unordered_map<std::string, std::string>;

inline FastaRef loadFasta(const std::string& path) {
    FastaRef ref;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "[BroCOLI] cannot open FASTA: " << path << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cout << "[BroCOLI] reading FASTA " << path << "\n";

    std::string line, name, seq;
    auto flush = [&] {
        if (!name.empty()) ref.emplace(name, std::move(seq));
        seq.clear();
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            size_t e = line.find_first_of(" \t");
            name = line.substr(1, e == std::string::npos ? std::string::npos : e - 1);
        } else {
            seq += line;
        }
    }
    flush();
    std::cout << "[BroCOLI] FASTA: " << ref.size() << " sequences\n";
    return ref;
}

// Canonical donor/acceptor check. Returns 1 (+ strand), 2 (- strand) or 0.
inline int spliceSignal(const std::string& seq, const Interval& sj, int strand) {
    if (sj[0] < 1 || sj[1] < 2) return 0;
    const size_t p1 = static_cast<size_t>(sj[0]) - 1;
    const size_t p2 = static_cast<size_t>(sj[1]) - 2;
    if (p1 + 2 > seq.size() || p2 + 2 > seq.size()) return 0;

    char d0 = std::toupper(static_cast<unsigned char>(seq[p1]));
    char d1 = std::toupper(static_cast<unsigned char>(seq[p1 + 1]));
    char a0 = std::toupper(static_cast<unsigned char>(seq[p2]));
    char a1 = std::toupper(static_cast<unsigned char>(seq[p2 + 1]));

    auto is = [](char x, char y, char c1, char c2) { return x == c1 && y == c2; };

    if (strand == 1) {
        if ((is(d0, d1, 'G', 'T') && is(a0, a1, 'A', 'G')) ||
            (is(d0, d1, 'G', 'C') && is(a0, a1, 'A', 'G')) ||
            (is(d0, d1, 'A', 'T') && is(a0, a1, 'A', 'C'))) return 1;
        return 0;
    }
    if ((is(d0, d1, 'C', 'T') && is(a0, a1, 'A', 'C')) ||
        (is(d0, d1, 'C', 'T') && is(a0, a1, 'G', 'C')) ||
        (is(d0, d1, 'G', 'T') && is(a0, a1, 'A', 'T'))) return 2;
    return 0;
}

// ------------------------------------------------------------ ref mapping --
// Different samples may order their @SQ lines differently, so every tid is
// translated into one global reference id before anything is written out.
class RefTable {
public:
    int add(const char* name, hts_pos_t len) {
        auto it = index_.find(name);
        if (it != index_.end()) return it->second;
        const int id = static_cast<int>(names_.size());
        names_.emplace_back(name);
        lens_.push_back(len);
        index_.emplace(names_.back(), id);
        return id;
    }
    int find(const std::string& name) const {
        auto it = index_.find(name);
        return it == index_.end() ? -1 : it->second;
    }
    const std::string& name(int id) const { return names_[id]; }
    hts_pos_t length(int id) const { return lens_[id]; }
    size_t size() const { return names_.size(); }

private:
    std::vector<std::string> names_;
    std::vector<hts_pos_t> lens_;
    std::unordered_map<std::string, int> index_;
};

// --------------------------------------------------------------- windows ---
struct Region { int tid; hts_pos_t beg, end; };   // BAM-local tid, 0-based half open
struct WindowTask {
    std::vector<Region> regions;                  // contiguous in genome order
    std::string out_path;
};

// Small contigs are chained into one task so that a reference with thousands of
// scaffolds does not produce thousands of tiny files.
inline std::vector<WindowTask> planWindows(sam_hdr_t* hdr, int window,
                                           const std::string& dir) {
    std::vector<WindowTask> tasks;
    WindowTask cur;
    hts_pos_t acc = 0;

    auto emit = [&] {
        if (cur.regions.empty()) return;
        char buf[64];
        std::snprintf(buf, sizeof(buf), "w%05zu.rec", tasks.size());
        cur.out_path = joinPath(dir, buf);
        tasks.push_back(std::move(cur));
        cur = WindowTask{};
        acc = 0;
    };

    const int nref = sam_hdr_nref(hdr);
    for (int tid = 0; tid < nref; ++tid) {
        hts_pos_t len = sam_hdr_tid2len(hdr, tid);
        if (len <= 0) len = INT32_MAX;
        for (hts_pos_t b = 0; b < len; b += window) {
            const hts_pos_t e = std::min<hts_pos_t>(b + window, len);
            cur.regions.push_back(Region{tid, b, e});
            acc += e - b;
            if (acc >= window) emit();
        }
    }
    emit();
    return tasks;
}

// ------------------------------------------------------- barcode / umi -----
inline std::string normaliseTag(const std::string& t) {
    return t.size() >= 2 ? t.substr(0, 2) : t;
}

// Flexiplex-style names look like "BARCODE-1_UMI#original_name".
inline bool parseFlexiplexName(const std::string& qname,
                               std::string& bc, std::string& umi) {
    const size_t hash = qname.find('#');
    if (hash == std::string::npos) return false;
    const size_t us = qname.find('_');
    if (us == std::string::npos || us > hash) return false;

    bc = qname.substr(0, us);
    if (bc.size() >= 2 && bc.compare(bc.size() - 2, 2, "-1") == 0)
        bc.resize(bc.size() - 2);
    umi = qname.substr(us + 1, hash - us - 1);
    return !bc.empty() && !umi.empty();
}

// ----------------------------------------------------------- stage 1 core --
struct Stage1Result {
    std::vector<std::string> chunk_files;   // in genome order
    std::set<std::string> barcodes;         // sc only
    long long n_records = 0;
};

class BamSampleReader {
public:
    BamSampleReader(const std::string& path, const RefTable& refs) : path_(path) {
        fp_ = sam_open(path.c_str(), "r");
        if (!fp_) throw std::runtime_error("cannot open BAM: " + path);
        hdr_ = sam_hdr_read(fp_);
        if (!hdr_) throw std::runtime_error("cannot read BAM header: " + path);
        idx_ = sam_index_load(fp_, path.c_str());
        if (!idx_) throw std::runtime_error("missing BAM index for " + path +
                                            " (run `samtools index`)");
        // local tid -> global ref id
        tid2gid_.resize(sam_hdr_nref(hdr_), -1);
        for (int i = 0; i < sam_hdr_nref(hdr_); ++i)
            tid2gid_[i] = refs.find(sam_hdr_tid2name(hdr_, i));
        rec_ = bam_init1();
    }
    ~BamSampleReader() {
        if (rec_) bam_destroy1(rec_);
        if (idx_) hts_idx_destroy(idx_);
        if (hdr_) sam_hdr_destroy(hdr_);
        if (fp_) hts_close(fp_);
    }
    BamSampleReader(const BamSampleReader&) = delete;
    BamSampleReader& operator=(const BamSampleReader&) = delete;

    sam_hdr_t* header() const { return hdr_; }
    int globalTid(int tid) const {
        return (tid >= 0 && tid < static_cast<int>(tid2gid_.size())) ? tid2gid_[tid] : -1;
    }

    hts_itr_t* query(const Region& r) { return sam_itr_queryi(idx_, r.tid, r.beg, r.end); }
    int next(hts_itr_t* it) { return sam_itr_next(fp_, it, rec_); }
    bam1_t* rec() const { return rec_; }

private:
    std::string path_;
    samFile* fp_ = nullptr;
    sam_hdr_t* hdr_ = nullptr;
    hts_idx_t* idx_ = nullptr;
    bam1_t* rec_ = nullptr;
    std::vector<int> tid2gid_;
};

// Converts one alignment. Returns false if the record must be skipped.
inline bool convertAlignment(const bam1_t* b, const Options& opt,
                             const std::string* chrom_seq,
                             const char* bc_tag, const char* umi_tag,
                             int global_tid, int sample_idx, ReadRec& out) {
    const bam1_core_t& c = b->core;

    if (c.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        ++g_stats.drop_flag; return false;
    }
    if (c.qual < opt.mapq) { ++g_stats.drop_mapq; return false; }
    if (c.n_cigar == 0 || global_tid < 0) { ++g_stats.drop_cigar; return false; }

    // ---- CIGAR -> aligned blocks (1-based start, end exclusive) ----
    const uint32_t* cig = bam_get_cigar(b);
    int pos = static_cast<int>(c.pos) + 1;
    int qlen = 0;
    IntervalVec blocks;
    blocks.reserve(c.n_cigar / 2 + 1);

    for (uint32_t i = 0; i < c.n_cigar; ++i) {
        const int op = bam_cigar_op(cig[i]);
        const int n  = static_cast<int>(bam_cigar_oplen(cig[i]));
        switch (op) {
            case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                blocks.push_back(Interval{pos, pos + n});
                pos += n; qlen += n; break;
            case BAM_CDEL: case BAM_CREF_SKIP:
                pos += n; break;
            case BAM_CINS: case BAM_CSOFT_CLIP:
                qlen += n; break;
            default: break;   // H, P
        }
    }
    if (blocks.empty()) { ++g_stats.drop_cigar; return false; }

    // Inherited consistency check: the CIGAR must account for the whole query.
    if (c.l_qseq > 0 && qlen != c.l_qseq) { ++g_stats.drop_cigar; return false; }

    // ---- barcode / UMI ----
    out.bc.clear();
    out.umi.clear();
    if (opt.mode == Mode::SC) {
        if (const uint8_t* p = bam_aux_get(b, bc_tag))
            if (const char* s = bam_aux2Z(p)) out.bc = s;
        if (const uint8_t* p = bam_aux_get(b, umi_tag))
            if (const char* s = bam_aux2Z(p)) out.umi = s;
        if (out.bc.empty() || out.umi.empty())
            parseFlexiplexName(bam_get_qname(b), out.bc, out.umi);
        if (out.bc.empty() || out.umi.empty()) { ++g_stats.drop_nobc; return false; }
        // 10x style "AACGTT-1" -> "AACGTT"
        if (out.bc.size() > 2 && out.bc.compare(out.bc.size() - 2, 2, "-1") == 0)
            out.bc.resize(out.bc.size() - 2);
    }

    out.name.assign(bam_get_qname(b));
    out.tid    = global_tid;
    out.file   = sample_idx;
    out.len    = c.l_qseq > 0 ? c.l_qseq : qlen;
    out.strand = (c.flag & BAM_FREVERSE) ? 0 : 1;

    // ---- splice junctions ----
    out.sjs.clear();
    out.sigs.clear();
    for (size_t i = 0; i + 1 < blocks.size(); ++i) {
        const int former = blocks[i][1];
        const int later  = blocks[i + 1][0];
        if (later - former <= opt.sj_distance) continue;

        out.sjs.push_back(Interval{former, later - 1});

        // Anchors shorter than 10 bp are not trusted (inherited rule).
        const int lanchor = (former - 1) - blocks[i][0];
        const int ranchor = blocks[i + 1][1] - later;
        int sig = 0;
        if (lanchor > 10 && ranchor > 10 && chrom_seq)
            sig = spliceSignal(*chrom_seq, out.sjs.back(), out.strand);
        out.sigs.push_back(sig);
    }

    // Note: 1.x reported `back()[1]` (one past the last aligned base) for
    // spliced reads and `back()[1] - 1` for single-exon reads. The inconsistency
    // leaked into novel isoform GTF coordinates, so both now use the last
    // aligned base.
    out.beg = blocks.front()[0];
    out.end = blocks.back()[1] - 1;
    ++g_stats.reads_kept;
    return true;
}

// Runs stage 1 for one sample. Returns the ordered list of chunk files.
inline Stage1Result runStage1Sample(const std::string& bam_path, int sample_idx,
                                    const Options& opt, const RefTable& refs,
                                    const FastaRef& fasta, const std::string& tmpdir,
                                    ThreadPool& pool) {
    Stage1Result result;

    // The window plan only needs the header.
    std::vector<WindowTask> tasks;
    {
        BamSampleReader probe(bam_path, refs);
        tasks = planWindows(probe.header(), opt.window, tmpdir);
    }

    const std::string bcTag  = normaliseTag(opt.bc_tag);
    const std::string umiTag = normaliseTag(opt.umi_tag);

    std::mutex mu;
    std::vector<char> nonempty(tasks.size(), 0);

    TaskQueue queue(pool);
    for (size_t ti = 0; ti < tasks.size(); ++ti) {
        queue.submit([&, ti] {
            const WindowTask& task = tasks[ti];
            BamSampleReader reader(bam_path, refs);
            RecordWriter writer(task.out_path);

            ReadRec rec;
            std::set<std::string> localBc;
            long long n = 0;
            int cachedTid = -1;
            const std::string* cachedSeq = nullptr;

            for (const Region& r : task.regions) {
                hts_itr_t* it = reader.query(r);
                if (!it) continue;
                while (reader.next(it) >= 0) {
                    const bam1_t* b = reader.rec();
                    ++g_stats.reads_in;
                    if (b->core.pos < r.beg) continue;      // owned by an earlier window

                    const int gid = reader.globalTid(b->core.tid);
                    if (gid != cachedTid) {
                        cachedTid = gid;
                        cachedSeq = nullptr;
                        if (gid >= 0) {
                            auto f = fasta.find(refs.name(gid));
                            if (f != fasta.end()) cachedSeq = &f->second;
                        }
                    }
                    if (!convertAlignment(b, opt, cachedSeq, bcTag.c_str(), umiTag.c_str(),
                                          gid, sample_idx, rec))
                        continue;

                    if (opt.mode == Mode::SC) localBc.insert(rec.bc);
                    writer.append(rec);
                    ++n;
                }
                hts_itr_destroy(it);
            }
            writer.close();

            std::lock_guard<std::mutex> lk(mu);
            nonempty[ti] = n > 0;
            result.n_records += n;
            if (opt.mode == Mode::SC)
                result.barcodes.insert(localBc.begin(), localBc.end());
        });
    }
    queue.wait();

    for (size_t ti = 0; ti < tasks.size(); ++ti) {
        if (nonempty[ti]) result.chunk_files.push_back(tasks[ti].out_path);
        else std::remove(tasks[ti].out_path.c_str());
    }
    return result;
}

}  // namespace brocoli
