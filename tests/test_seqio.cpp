// preBroCOLI 2.0 - sequence I/O tests. Needs zlib only (no edlib, no WFA2).
#include "seqio.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace prebrocoli;

namespace {

int g_fail = 0;

void check(bool ok, const std::string& what) {
    if (!ok) { std::cout << "FAIL: " << what << "\n"; ++g_fail; }
}

std::string tmp(const std::string& name) { return "/tmp/prebrocoli_test_" + name; }

void writeFile(const std::string& path, const std::string& content) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    std::fwrite(content.data(), 1, content.size(), f);
    std::fclose(f);
}

void gzipFile(const std::string& src, const std::string& dst) {
    const std::string cmd = "gzip -c " + src + " > " + dst;
    if (std::system(cmd.c_str()) != 0) { std::cout << "FAIL: gzip failed\n"; ++g_fail; }
}

std::vector<SeqRecord> readAll(const std::string& path) {
    SeqReader r(path);
    SeqRecord rec;
    std::vector<SeqRecord> out;
    while (r.next(rec)) out.push_back(rec);
    return out;
}

const char* kFastq =
    "@read1 runid=abc len=8\n"
    "ACGTACGT\n"
    "+\n"
    "IIIIIIII\n"
    "@read2\n"
    "TTTTAAAA\n"
    "+read2\n"
    "01234567\n";

void testPlainFastq() {
    const std::string p = tmp("a.fastq");
    writeFile(p, kFastq);
    const auto recs = readAll(p);
    check(recs.size() == 2, "plain fastq: two records");
    check(recs[0].id == "read1", "plain fastq: id stops at first whitespace");
    check(recs[0].seq == "ACGTACGT", "plain fastq: sequence");
    check(recs[0].qual == "IIIIIIII", "plain fastq: quality");
    check(recs[1].id == "read2", "plain fastq: second id");
    check(recs[1].qual == "01234567", "plain fastq: second quality");
    SeqReader r(p);
    check(r.isFastq(), "plain fastq: detected as fastq");
}

void testGzMatchesPlain() {
    const std::string p = tmp("a.fastq"), g = tmp("a.fastq.gz");
    gzipFile(p, g);
    const auto a = readAll(p);
    const auto b = readAll(g);
    check(a.size() == b.size(), "gz: same record count as plain");
    bool same = a.size() == b.size();
    for (size_t i = 0; same && i < a.size(); ++i)
        same = a[i].id == b[i].id && a[i].seq == b[i].seq && a[i].qual == b[i].qual;
    check(same, "gz: records identical to plain, same code path");
}

void testCrlf() {
    const std::string p = tmp("crlf.fastq");
    writeFile(p, "@r1\r\nACGT\r\n+\r\nIIII\r\n");
    const auto recs = readAll(p);
    check(recs.size() == 1, "crlf: one record");
    check(recs[0].seq == "ACGT", "crlf: carriage return stripped from sequence");
    check(recs[0].qual == "IIII", "crlf: carriage return stripped from quality");
}

void testMultiLineFasta() {
    const std::string p = tmp("a.fasta");
    writeFile(p, ">r1 desc\nACGT\nACGT\nAC\n>r2\nTTTT\n");
    const auto recs = readAll(p);
    check(recs.size() == 2, "fasta: two records");
    check(recs[0].seq == "ACGTACGTAC", "fasta: wrapped sequence joined");
    check(recs[0].qual.empty(), "fasta: no quality");
    SeqReader r(p);
    check(!r.isFastq(), "fasta: detected as fasta");
}

void testEmpty() {
    const std::string p = tmp("empty.fastq");
    writeFile(p, "");
    const auto recs = readAll(p);
    check(recs.empty(), "empty file: no records, no throw");
}

void testTruncated() {
    const std::string p = tmp("trunc.fastq");
    writeFile(p, "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n");
    bool threw = false;
    try { readAll(p); } catch (const std::exception&) { threw = true; }
    check(threw, "truncated record: throws instead of silently producing garbage");
}

void testLengthMismatch() {
    const std::string p = tmp("mismatch.fastq");
    writeFile(p, "@r1\nACGTA\n+\nIIII\n");
    bool threw = false;
    try { readAll(p); } catch (const std::exception&) { threw = true; }
    check(threw, "seq/qual length mismatch: throws");
}

void testMissingFile() {
    bool threw = false;
    try { SeqReader r("/tmp/prebrocoli_test_does_not_exist.fastq"); }
    catch (const std::exception&) { threw = true; }
    check(threw, "missing file: throws");
}

void testRoundTrip() {
    const std::string p = tmp("out.fastq");
    {
        LineSink sink(p);
        for (int i = 0; i < 5000; ++i)
            sink.record("read" + std::to_string(i) + "\tCB:Z:ACGT\tUB:Z:TTTT",
                        std::string(120, 'A'), std::string(120, 'I'), true);
    }
    const auto recs = readAll(p);
    check(recs.size() == 5000, "writer: all records survive buffering");
    check(recs[4999].id == "read4999", "writer: id token stops at tab");
    check(recs[0].seq.size() == 120 && recs[0].qual.size() == 120,
          "writer: sequence and quality round-trip");
}

void testLargeBuffered() {
    // Crosses the reader's 1 MB refill boundary many times.
    const std::string p = tmp("big.fastq");
    std::string content;
    content.reserve(6u << 20);
    for (int i = 0; i < 20000; ++i) {
        const std::string s(200, "ACGT"[i % 4]);
        content += "@r" + std::to_string(i) + "\n" + s + "\n+\n" + std::string(200, 'I') + "\n";
    }
    writeFile(p, content);
    const std::string g = tmp("big.fastq.gz");
    gzipFile(p, g);
    const auto a = readAll(p);
    const auto b = readAll(g);
    check(a.size() == 20000, "large input: record count");
    check(a.size() == b.size() && a[19999].seq == b[19999].seq,
          "large input: gz and plain agree across buffer refills");
}

}  // namespace

int main() {
    testPlainFastq();
    testGzMatchesPlain();
    testCrlf();
    testMultiLineFasta();
    testEmpty();
    testTruncated();
    testLengthMismatch();
    testMissingFile();
    testRoundTrip();
    testLargeBuffered();

    if (g_fail == 0) std::cout << "ALL SEQIO TESTS PASSED (0 failures)\n";
    else std::cout << g_fail << " FAILURES\n";
    return g_fail == 0 ? 0 : 1;
}
