// Intermediate record format.
//
// BroCOLI 1.x wrote a tab-separated text file between stage 1 and stage 2 and
// re-parsed it with istringstream + stoi (three allocations and a dozen
// integer parses per read). We keep the same two-stage design - it is what
// bounds memory and what lets several samples be merged - but the payload is
// now a packed binary record, so stage 2 is a memcpy instead of a parse.
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#include "common.hpp"

namespace brocoli {

#pragma pack(push, 1)
struct RecHeader {
    int32_t  tid;        // reference id (index into the BAM header)
    int32_t  beg;        // read start, 1-based
    int32_t  end;        // read end (see note in README about the inherited off-by-one)
    int32_t  len;        // query length
    uint16_t file;       // sample index
    uint16_t n_sj;
    uint16_t name_len;
    uint8_t  bc_len;
    uint8_t  umi_len;
    uint8_t  strand;     // 1 = forward, 0 = reverse
    uint8_t  pad{0};
};
#pragma pack(pop)

struct RecSj {
    int32_t beg;
    int32_t end;
    int32_t sig;         // 0 = unsupported, 1 = canonical +, 2 = canonical -
};

// Fully materialised record.
struct ReadRec {
    std::string name, bc, umi;
    int tid = -1, beg = 0, end = 0, len = 0;
    int file = 0, strand = 1;
    IntervalVec sjs;
    std::vector<int> sigs;
};

// ------------------------------------------------------------------ write --
class RecordWriter {
public:
    explicit RecordWriter(const std::string& path) : path_(path) {
        fp_ = std::fopen(path.c_str(), "wb");
        if (!fp_) throw std::runtime_error("cannot open for writing: " + path);
        buf_.reserve(kFlush + (1 << 16));
    }
    ~RecordWriter() { close(); }

    RecordWriter(const RecordWriter&) = delete;
    RecordWriter& operator=(const RecordWriter&) = delete;

    void append(const ReadRec& r) {
        RecHeader h{};
        h.tid = r.tid; h.beg = r.beg; h.end = r.end; h.len = r.len;
        h.file = static_cast<uint16_t>(r.file);
        h.n_sj = static_cast<uint16_t>(r.sjs.size());
        h.name_len = static_cast<uint16_t>(r.name.size());
        h.bc_len   = static_cast<uint8_t>(r.bc.size());
        h.umi_len  = static_cast<uint8_t>(r.umi.size());
        h.strand   = static_cast<uint8_t>(r.strand);

        raw(&h, sizeof(h));
        raw(r.name.data(), r.name.size());
        raw(r.bc.data(), r.bc.size());
        raw(r.umi.data(), r.umi.size());
        for (size_t i = 0; i < r.sjs.size(); ++i) {
            RecSj s{r.sjs[i][0], r.sjs[i][1], i < r.sigs.size() ? r.sigs[i] : 0};
            raw(&s, sizeof(s));
        }
        if (buf_.size() >= kFlush) flush();
    }

    // Copy an already-encoded record straight through (used by the merger).
    void appendRaw(const char* p, size_t n) {
        raw(p, n);
        if (buf_.size() >= kFlush) flush();
    }

    void flush() {
        if (buf_.empty()) return;
        if (std::fwrite(buf_.data(), 1, buf_.size(), fp_) != buf_.size())
            throw std::runtime_error("short write on " + path_);
        written_ += static_cast<long long>(buf_.size());
        buf_.clear();
    }

    // Current offset, tracked by hand. Calling ftell() here instead would force
    // a flush at every group boundary and defeat the write buffer entirely.
    long long tell() const { return written_ + static_cast<long long>(buf_.size()); }

    void close() {
        if (!fp_) return;
        flush();
        std::fclose(fp_);
        fp_ = nullptr;
    }

private:
    void raw(const void* p, size_t n) {
        if (n) buf_.append(static_cast<const char*>(p), n);
    }
    static constexpr size_t kFlush = 8u << 20;
    long long written_ = 0;
    std::string path_, buf_;
    std::FILE* fp_ = nullptr;
};

// ------------------------------------------------------------------- read --
class RecordReader {
public:
    explicit RecordReader(const std::string& path) : path_(path) {
        fp_ = std::fopen(path.c_str(), "rb");
        if (!fp_) throw std::runtime_error("cannot open for reading: " + path);
        std::setvbuf(fp_, nullptr, _IOFBF, 4u << 20);
    }
    ~RecordReader() { if (fp_) std::fclose(fp_); }

    RecordReader(const RecordReader&) = delete;
    RecordReader& operator=(const RecordReader&) = delete;

    void seek(long long off) { std::fseek(fp_, off, SEEK_SET); }

    // Reads one record. `blob`, when non-null, receives the exact bytes so the
    // merger can forward them without re-encoding.
    bool next(ReadRec& r, std::string* blob = nullptr) {
        RecHeader h;
        if (std::fread(&h, 1, sizeof(h), fp_) != sizeof(h)) return false;

        const size_t payload = size_t(h.name_len) + h.bc_len + h.umi_len +
                               size_t(h.n_sj) * sizeof(RecSj);
        scratch_.resize(payload);
        if (payload && std::fread(scratch_.data(), 1, payload, fp_) != payload)
            throw std::runtime_error("truncated record in " + path_);

        r.tid = h.tid; r.beg = h.beg; r.end = h.end; r.len = h.len;
        r.file = h.file; r.strand = h.strand;

        const char* p = scratch_.data();
        r.name.assign(p, h.name_len); p += h.name_len;
        r.bc.assign(p, h.bc_len);     p += h.bc_len;
        r.umi.assign(p, h.umi_len);   p += h.umi_len;

        r.sjs.clear(); r.sigs.clear();
        r.sjs.reserve(h.n_sj); r.sigs.reserve(h.n_sj);
        for (int i = 0; i < h.n_sj; ++i) {
            RecSj s;
            std::memcpy(&s, p, sizeof(s));
            p += sizeof(s);
            r.sjs.push_back(Interval{s.beg, s.end});
            r.sigs.push_back(s.sig);
        }

        if (blob) {
            blob->assign(reinterpret_cast<const char*>(&h), sizeof(h));
            blob->append(scratch_.data(), payload);
        }
        return true;
    }

private:
    std::string path_;
    std::vector<char> scratch_;
    std::FILE* fp_ = nullptr;
};

// ------------------------------------------------------------ group index --
struct GroupIndexEntry {
    long long offset = 0;    // byte offset into the merged record file
    int  n_reads = 0;
    int  tid = -1;
    int  low = 0, high = 0;
    int  gid = 0;            // 1-based, used when naming novel isoforms
};

}  // namespace brocoli
