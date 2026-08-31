// Shared timestamped logging for BroCOLI 2.0 and preBroCOLI 2.0.
//
// Extracted verbatim from brocoli/src/main.cpp so both tools produce the same
// log format. Header-only, depends on nothing but the standard library.
#ifndef BROCOLI_LOG_HPP
#define BROCOLI_LOG_HPP

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <mutex>
#include <streambuf>
#include <string>

namespace brolog {

// --------------------------------------------------------------- logging --
//
// Two problems are fixed here at once.
//
// 1. Interleaving. 1.x (and 2.0 up to now) called sync_with_stdio(false) and
//    then mixed three sinks: std::cout (iostream's own buffer, fully buffered
//    when redirected to a file), std::printf (C stdout buffer) and
//    std::fprintf(stderr) (unbuffered). Redirecting to a log file produced
//    output in an order that had nothing to do with the order things were
//    printed - the banner could show up *after* the per-cluster progress
//    lines, and the region summary before it. There is now exactly one path
//    to the terminal: everything goes through std::cout, and this buffer
//    flushes the destination at every newline, so the file order is the
//    print order. sync_with_stdio is deliberately left at its default (on)
//    so that any stray printf still lands in the same stream in the right
//    place. Nothing here is on a hot path - it is a few hundred lines total.
//
// 2. Timestamps. Every line gets a "2026-08-30 18:54:08,392 - INFO - " prefix.
//    Doing it in the streambuf rather than at each call site means log lines
//    printed from the other headers (gtf.hpp, bam_stage1.hpp, group.hpp) are
//    stamped too, without touching those files.
//
// The prefix is inserted at the start of each line, so partial writes are
// fine; the mutex is shared by the cout and cerr filters so that lines from
// worker threads cannot be spliced into each other.
class LogBuf : public std::streambuf {
public:
    LogBuf(std::streambuf* dest, const char* level, std::mutex& mtx)
        : dest_(dest), level_(level), mtx_(mtx) {}

protected:
    int overflow(int ch) override {
        if (ch == traits_type::eof()) return traits_type::not_eof(ch);
        const char c = static_cast<char>(ch);
        std::lock_guard<std::mutex> lock(mtx_);
        return emit(&c, 1) ? ch : traits_type::eof();
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        std::lock_guard<std::mutex> lock(mtx_);
        return emit(s, static_cast<size_t>(n)) ? n : 0;
    }

    int sync() override {
        std::lock_guard<std::mutex> lock(mtx_);
        return dest_->pubsync();
    }

private:
    static std::string stamp() {
        using namespace std::chrono;
        const auto now = system_clock::now();
        const std::time_t secs = system_clock::to_time_t(now);
        const auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
        std::tm tmv{};
        localtime_r(&secs, &tmv);
        char date[24];
        std::strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S", &tmv);
        char out[40];
        std::snprintf(out, sizeof(out), "%s,%03d", date,
                      static_cast<int>(ms.count()));
        return out;
    }

    bool emit(const char* s, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            if (line_start_) {
                const std::string p = stamp() + " - " + level_ + " - ";
                if (dest_->sputn(p.data(), static_cast<std::streamsize>(p.size()))
                    != static_cast<std::streamsize>(p.size())) return false;
                line_start_ = false;
            }
            if (dest_->sputc(s[i]) == traits_type::eof()) return false;
            if (s[i] == '\n') {
                line_start_ = true;
                dest_->pubsync();   // keeps file order == print order
            }
        }
        return true;
    }

    std::streambuf* dest_;
    const char* level_;
    std::mutex& mtx_;
    bool line_start_ = true;
};

// Installs the filters on cout/cerr and puts the original buffers back on the
// way out, so the streams are never left pointing at a dead object.
class LogRedirect {
public:
    LogRedirect()
        : out_(std::cout.rdbuf(), "INFO", mtx_),
          err_(std::cerr.rdbuf(), "ERROR", mtx_),
          old_out_(std::cout.rdbuf(&out_)),
          old_err_(std::cerr.rdbuf(&err_)) {}

    ~LogRedirect() {
        std::cout.flush();
        std::cerr.flush();
        std::cout.rdbuf(old_out_);
        std::cerr.rdbuf(old_err_);
    }

    LogRedirect(const LogRedirect&) = delete;
    LogRedirect& operator=(const LogRedirect&) = delete;

private:
    std::mutex mtx_;
    LogBuf out_, err_;
    std::streambuf* old_out_;
    std::streambuf* old_err_;
};

// Formats one complete line and writes it with a single stream operation, so
// that concurrent writers cannot interleave mid-line.
void logLine(const std::string& msg) { std::cout << (msg + "\n"); }

std::string fmt(const char* pattern, ...) __attribute__((format(printf, 1, 2)));
std::string fmt(const char* pattern, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, pattern);
    std::vsnprintf(buf, sizeof(buf), pattern, ap);
    va_end(ap);
    return buf;
}

}  // namespace brolog

#endif  // BROCOLI_LOG_HPP
