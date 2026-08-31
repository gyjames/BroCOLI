// preBroCOLI 2.0 - sequence input / output.
//
// Input:  SeqReader reads plain FASTQ/FASTA, gzip-compressed FASTQ/FASTA and
//         stdin through one interface. zlib's gzread passes uncompressed data
//         through transparently, so a single code path covers .fastq and
//         .fastq.gz and no filename sniffing is needed.
//
// Output: LineSink is a buffered FILE* writer for plain text. preBroCOLI only
//         ever writes uncompressed .fastq / .tsv - see README for why.
#ifndef PREBROCOLI_SEQIO_HPP
#define PREBROCOLI_SEQIO_HPP

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <zlib.h>

namespace prebrocoli {

struct SeqRecord {
    std::string id;    // first whitespace-delimited token, without @ or >
    std::string seq;
    std::string qual;  // empty for FASTA
};

class SeqReader {
public:
    // path may be "-" or empty for stdin.
    explicit SeqReader(const std::string& path) {
        if (path.empty() || path == "-") {
            fp_ = gzdopen(fileno(stdin), "rb");
            name_ = "<stdin>";
        } else {
            fp_ = gzopen(path.c_str(), "rb");
            name_ = path;
        }
        if (!fp_) throw std::runtime_error("cannot open " + name_);
        gzbuffer(fp_, 1u << 20);          // 1 MB inflate buffer
        buf_.resize(1u << 20);
    }

    ~SeqReader() { if (fp_) gzclose(fp_); }

    SeqReader(const SeqReader&) = delete;
    SeqReader& operator=(const SeqReader&) = delete;

    // Priming reads the first header line, so this is valid before next().
    bool isFastq() { prime(); return fastq_; }
    const std::string& name() const { return name_; }

    // Fills rec and returns true, or returns false at end of input.
    bool next(SeqRecord& rec) {
        prime();
        if (header_.empty()) return false;

        rec.id = firstToken(header_, 1);
        rec.seq.clear();
        rec.qual.clear();

        if (fastq_) {
            if (!readLine(rec.seq)) throw std::runtime_error("truncated record: " + rec.id);
            std::string plus;
            if (!readLine(plus)) throw std::runtime_error("truncated record: " + rec.id);
            if (!readLine(rec.qual)) throw std::runtime_error("truncated record: " + rec.id);
            if (rec.qual.size() != rec.seq.size())
                throw std::runtime_error("seq/qual length mismatch in read " + rec.id);
            if (!readLine(header_)) header_.clear();
        } else {
            // FASTA: sequence may be wrapped over any number of lines.
            std::string chunk;
            header_.clear();
            while (readLine(chunk)) {
                if (!chunk.empty() && chunk[0] == '>') { header_ = chunk; break; }
                rec.seq += chunk;
            }
        }
        return true;
    }

private:
    // Finds the first header line and decides FASTA vs FASTQ. Idempotent.
    void prime() {
        if (primed_) return;
        primed_ = true;
        while (readLine(header_) && header_.empty()) {}
        if (header_.empty()) return;                       // empty input
        if (header_[0] == '@')      fastq_ = true;
        else if (header_[0] == '>') fastq_ = false;
        else throw std::runtime_error("not FASTA/FASTQ: " + name_);
    }

    static std::string firstToken(const std::string& line, size_t from) {
        size_t e = from;
        while (e < line.size() && line[e] != ' ' && line[e] != '\t') ++e;
        return line.substr(from, e - from);
    }

    // Buffered line reader. gzgets() would copy twice and caps at the caller's
    // buffer size; this refills a 1 MB block and hands out lines from it.
    bool readLine(std::string& out) {
        out.clear();
        for (;;) {
            if (pos_ == end_) {
                const int n = gzread(fp_, &buf_[0], static_cast<unsigned>(buf_.size()));
                if (n < 0) {
                    int err = 0;
                    const char* msg = gzerror(fp_, &err);
                    throw std::runtime_error("read error on " + name_ + ": " +
                                             (msg ? msg : "unknown"));
                }
                if (n == 0) return !out.empty();
                pos_ = 0;
                end_ = static_cast<size_t>(n);
            }
            const char* nl = static_cast<const char*>(
                std::memchr(buf_.data() + pos_, '\n', end_ - pos_));
            if (nl) {
                const size_t len = static_cast<size_t>(nl - (buf_.data() + pos_));
                out.append(buf_.data() + pos_, len);
                pos_ += len + 1;
                if (!out.empty() && out.back() == '\r') out.pop_back();  // CRLF
                return true;
            }
            out.append(buf_.data() + pos_, end_ - pos_);
            pos_ = end_;
        }
    }

    gzFile fp_ = nullptr;
    std::string name_;
    std::string buf_;
    std::string header_;
    size_t pos_ = 0, end_ = 0;
    bool primed_ = false;
    bool fastq_ = true;
};

// ------------------------------------------------------------------ output --
//
// Plain text only, on purpose. The intermediate FASTQ goes straight into
// minimap2 and is usually deleted afterwards; gzip costs about 3x the write
// time for roughly half the bytes, which is a bad trade for a scratch file.
class LineSink {
public:
    LineSink() = default;

    explicit LineSink(const std::string& path) { open(path); }

    ~LineSink() { close(); }

    LineSink(const LineSink&) = delete;
    LineSink& operator=(const LineSink&) = delete;

    void open(const std::string& path) {
        close();
        path_ = path;
        fp_ = std::fopen(path.c_str(), "wb");
        if (!fp_) throw std::runtime_error("cannot write " + path);
        buf_.reserve(kFlush + (1u << 16));
    }

    bool isOpen() const { return fp_ != nullptr; }
    const std::string& path() const { return path_; }

    void append(const std::string& s) {
        buf_ += s;
        if (buf_.size() >= kFlush) flush();
    }

    // One FASTQ / FASTA record.
    void record(const std::string& id, const std::string& seq,
                const std::string& qual, bool fastq) {
        buf_ += fastq ? '@' : '>';
        buf_ += id;
        buf_ += '\n';
        buf_ += seq;
        buf_ += '\n';
        if (fastq) {
            buf_ += '+';
            buf_ += id;
            buf_ += '\n';
            buf_ += qual;
            buf_ += '\n';
        }
        if (buf_.size() >= kFlush) flush();
    }

    void flush() {
        if (buf_.empty() || !fp_) return;
        if (std::fwrite(buf_.data(), 1, buf_.size(), fp_) != buf_.size())
            throw std::runtime_error("short write on " + path_);
        buf_.clear();
    }

    void close() {
        if (!fp_) return;
        flush();
        std::fclose(fp_);
        fp_ = nullptr;
    }

private:
    static constexpr size_t kFlush = 4u << 20;
    std::FILE* fp_ = nullptr;
    std::string path_;
    std::string buf_;
};

}  // namespace prebrocoli

#endif  // PREBROCOLI_SEQIO_HPP
