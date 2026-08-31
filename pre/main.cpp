// preBroCOLI 2.0 - barcode / UMI extraction from long-read FASTQ.
//
// Input  : .fastq or .fastq.gz (or "-" for stdin), auto-detected.
// Output : plain files under --outdir. Nothing goes to stdout any more, so
//          no shell redirection is needed and the log cannot end up mixed
//          into the data.
#include <getopt.h>

#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "log.hpp"

#include "common.hpp"
#include "demux.hpp"
#include "seqio.hpp"
#include "whitelist.hpp"

using namespace prebrocoli;
using brolog::fmt;
using brolog::logLine;
using brolog::LogRedirect;

namespace {

const char* kBanner = R"(
     ____               ____   ___  _      ___
    | __ )  _ __ ___   / ___| / _ \| |    |_ _|
    |  _ \ | '__/ _ \ | |    | | | | |     | |
    | |_) || | | (_) || |___ | |_| | |___  | |
    |____/ |_|  \___/  \____| \___/|_____||___|
)";

struct option kLongOpts[] = {
    {"chemistry",      required_argument, nullptr, 'q'},
    {"reads",          required_argument, nullptr, 'I'},
    {"whitelist",      required_argument, nullptr, 'w'},
    {"barcodeX",       required_argument, nullptr, 'x'},
    {"barcodeY",       required_argument, nullptr, 'y'},
    {"barcodeZ",       required_argument, nullptr, 'z'},
    {"outdir",         required_argument, nullptr, 'o'},
    {"prefix",         required_argument, nullptr, 'n'},
    {"thread",         required_argument, nullptr, 'p'},
    {"flank_editd",    required_argument, nullptr, 'f'},
    {"trim",           required_argument, nullptr, 'i'},
    {"split",          required_argument, nullptr, 's'},
    {"chimeric",       required_argument, nullptr, 'c'},
    {"batch",          required_argument, nullptr, 1001},
    {"keep_unmatched", no_argument,       nullptr, 1002},
    {"help",           no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

void printUsage(const char* prog) {
    std::cout
        << "Usage: " << prog << " -q <chemistry> -w <whitelist> -o <outdir> <reads.fastq[.gz]>\n\n"
        << "  reads may be .fastq, .fastq.gz, .fasta, .fasta.gz, or - for stdin.\n\n"
        << "Required:\n"
        << "  -q, --chemistry     magicseq | visium        (10x3v3 not implemented yet)\n"
        << "  -o, --outdir        output directory\n\n"
        << "Barcode source (one of):\n"
        << "  -w, --whitelist     whitelist file: barcode in column 1, optional UMI after it\n"
        << "  -x, --barcodeX      barcode list (magicseq: first segment)\n"
        << "  -y, --barcodeY      magicseq only, second segment\n"
        << "  -z, --barcodeZ      magicseq only, third segment\n\n"
        << "Common:\n"
        << "  -I, --reads         reads file (alternative to the positional argument)\n"
        << "  -n, --prefix        output filename prefix          (default: preBroCOLI)\n"
        << "  -p, --thread        worker threads                  (default: 1)\n"
        << "  -f, --flank_editd   max adapter edit distance       (magicseq 20, visium 8)\n"
        << "  -i, --trim          true/false: rewrite read ID and trim adapters (default: true)\n"
        << "  -s, --split         true/false: one file per barcode           (default: false)\n"
        << "  -c, --chimeric      true/false: mark chimeric reads with _C    (default: false)\n"
        << "      --keep_unmatched  also write reads with no barcode        (default: off)\n"
        << "      --batch         reads per thread per batch      (default: 2000)\n"
        << "  -h, --help\n\n"
        << "Outputs (under --outdir):\n"
        << "  <prefix>_matched.fastq          reads with a barcode\n"
        << "  <prefix>_unmatched.fastq        only with --keep_unmatched\n"
        << "  <prefix>_reads_barcodes.tsv     per-read barcode / UMI / edit distances\n"
        << "  <prefix>_barcode_counts.tsv     barcode -> read count\n"
        << "  <prefix>_summary.txt            same numbers as the terminal summary\n";
}

bool boolArg(const std::string& v) {
    std::string s = v;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "true" || s == "t" || s == "1") return true;
    if (s == "false" || s == "f" || s == "0") return false;
    std::cerr << "[preBroCOLI] expected true/false, got '" << v << "'\n";
    std::exit(1);
}

// ------------------------------------------------------------------ output --
class OutputWriter {
public:
    OutputWriter(const Options& opt, bool fastq) : opt_(opt), fastq_(fastq) {
        const std::string ext = fastq ? ".fastq" : ".fasta";
        matched_.open(joinPath(opt.outdir, opt.prefix + "_matched" + ext));
        if (opt.keep_unmatched)
            unmatched_.open(joinPath(opt.outdir, opt.prefix + "_unmatched" + ext));
        stats_.open(joinPath(opt.outdir, opt.prefix + "_reads_barcodes.tsv"));
        stats_.append("Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI\n");
    }

    void close() {
        matched_.close();
        unmatched_.close();
        stats_.close();
        for (auto& kv : split_) kv.second->close();
    }

    const std::map<std::string, long long>& barcodeCounts() const { return counts_; }

    void write(SearchResult& r) {
        for (const auto& b : r.vec_bc_for) counts_[b.barcode]++;
        for (const auto& b : r.vec_bc_rev) counts_[b.barcode]++;

        if (r.count == 0) {
            bump(g_stats.reads_unmatched);
            if (opt_.keep_unmatched && unmatched_.isOpen()) {
                unmatched_.record(r.read_id, r.line, r.qual_scores, fastq_);
                bump(g_stats.reads_written);
            }
            return;
        }

        bump(g_stats.reads_assigned);
        if (r.chimeric) bump(g_stats.reads_chimeric);
        if (r.count > 1) bump(g_stats.reads_multi_barcode);

        writeStats(r.read_id, r.vec_bc_for);
        writeStats(r.read_id, r.vec_bc_rev);

        emit(r.read_id + "_+", r.line, r.qual_scores, r.vec_bc_for,
             opt_.mark_chimeric && r.chimeric);

        if (opt_.trim_barcodes || r.vec_bc_for.empty()) {
            std::reverse(r.qual_scores.begin(), r.qual_scores.end());
            emit(r.read_id + "_-", r.rev_line, r.qual_scores, r.vec_bc_rev,
                 opt_.mark_chimeric && r.chimeric);
        }
    }

private:
    void writeStats(const std::string& read_id, const std::vector<Barcode>& v) {
        for (const auto& b : v)
            stats_.append(read_id + "\t" + b.barcode + "\t" +
                          std::to_string(b.flank_editd) + "\t" +
                          std::to_string(b.editd) + "\t" + b.umi + "\n");
    }

    LineSink& sinkFor(const std::string& barcode) {
        if (!opt_.split_by_barcode) return matched_;
        auto it = split_.find(barcode);
        if (it != split_.end()) return *it->second;
        // 1.x opened and closed a stream for every single read here.
        std::unique_ptr<LineSink> s(new LineSink(
            joinPath(opt_.outdir, opt_.prefix + "_" + barcode +
                                      (fastq_ ? ".fastq" : ".fasta"))));
        LineSink& ref = *s;
        split_.emplace(barcode, std::move(s));
        return ref;
    }

    // Port of 1.x print_read: trims everything up to the end of the adapter
    // and rewrites the read ID as barcode_umi#readid<i>of<n>[_C] plus the
    // CB:Z:/UB:Z: tags.
    void emit(const std::string& read_id, const std::string& read,
              const std::string& qual, const std::vector<Barcode>& vec_bc,
              bool chimeric) {
        const size_t vec_size = vec_bc.size();
        for (size_t b = 0; b < vec_size; ++b) {
            std::string suffix = std::to_string(b + 1) + "of" + std::to_string(vec_size);
            if (chimeric) suffix += "_C";

            const std::string& barcode = vec_bc[b].barcode;
            std::string new_read_id = barcode + "_" + vec_bc[b].umi + "#" + read_id +
                                      suffix + "\tCB:Z:" + barcode + "\tUB:Z:" + vec_bc[b].umi;

            int read_start = vec_bc[b].flank_end + 1;
            int read_length = static_cast<int>(read.length()) - read_start;
            for (size_t f = 0; f < vec_size; ++f) {
                const int t = vec_bc[f].flank_start - read_start;
                if (t >= 0 && t < read_length) read_length = t;
            }

            std::string qual_new, read_new;
            bool last = false;
            if (b == 0 && !opt_.trim_barcodes) {
                new_read_id = read_id;
                read_new = read;
                qual_new = qual;
                last = true;                 // 1.x forced the loop to exit here
            } else {
                if (read_start < 0 || read_length <= 0 ||
                    static_cast<size_t>(read_start + read_length) > read.length())
                    continue;
                if (!qual.empty()) {
                    if (static_cast<size_t>(read_start + read_length) > qual.length()) {
                        std::cerr << "[preBroCOLI] sequence and quality lengths differ for read "
                                  << read_id << ", skipping\n";
                        return;
                    }
                    qual_new = qual.substr(static_cast<size_t>(read_start),
                                           static_cast<size_t>(read_length));
                }
                read_new = read.substr(static_cast<size_t>(read_start),
                                       static_cast<size_t>(read_length));
            }

            if (read_new.empty()) continue;
            sinkFor(barcode).record(new_read_id, read_new, qual_new, fastq_);
            bump(g_stats.reads_written);
            if (last) break;
        }
    }

    const Options& opt_;
    bool fastq_;
    LineSink matched_, unmatched_, stats_;
    std::map<std::string, std::unique_ptr<LineSink>> split_;
    std::map<std::string, long long> counts_;
};

std::string summaryText(const Options& opt, double wall) {
    const long long in = g_stats.reads_in.load();
    auto pct = [&](long long v) {
        return in ? 100.0 * static_cast<double>(v) / static_cast<double>(in) : 0.0;
    };
    std::string s;
    s += "-----------------------------------------------------\n";
    s += "Summary:\n";
    s += fmt("Chemistry:             %s\n", opt.mode.c_str());
    s += fmt("Reads in:              %lld\n", in);
    s += fmt("Reads with barcode:    %lld (%.2f%%)\n",
             g_stats.reads_assigned.load(), pct(g_stats.reads_assigned.load()));
    s += fmt("Reads without barcode: %lld (%.2f%%)%s\n",
             g_stats.reads_unmatched.load(), pct(g_stats.reads_unmatched.load()),
             opt.keep_unmatched ? "" : "  [discarded, use --keep_unmatched to save]");
    s += fmt("Reads chimeric:        %lld\n", g_stats.reads_chimeric.load());
    s += fmt("Reads >1 barcode:      %lld\n", g_stats.reads_multi_barcode.load());
    s += fmt("Records written:       %lld\n", g_stats.reads_written.load());
    s += "-- where barcode searches failed (per search attempt, both strands) --\n";
    s += fmt("read_too_short:        %lld\n", g_stats.read_too_short.load());
    s += fmt("adapter_not_found:     %lld\n", g_stats.adapter_not_found.load());
    s += fmt("bc1_unresolved:        %lld\n", g_stats.bc1_unresolved.load());
    s += fmt("bc2_unresolved:        %lld\n", g_stats.bc2_unresolved.load());
    s += fmt("bc3_unresolved:        %lld\n", g_stats.bc3_unresolved.load());
    s += fmt("whitelist_no_match:    %lld\n", g_stats.whitelist_no_match.load());
    s += fmt("no_confident_call:     %lld\n", g_stats.search_no_call.load());
    s += fmt("barcode_editd_high:    %lld\n", g_stats.barcode_editd_high.load());
    s += "-- UMI --\n";
    s += fmt("umi_exact:             %lld\n", g_stats.umi_exact.load());
    s += fmt("umi_corrected:         %lld\n", g_stats.umi_corrected.load());
    s += fmt("umi_uncorrected:       %lld\n", g_stats.umi_uncorrected.load());
    s += "-- CPU time (summed over threads, not wall time) --\n";
    s += fmt("adapter search:        %.1f s\n", g_stats.t_adapter_us.load() / 1e6);
    s += fmt("barcode 1:             %.1f s\n", g_stats.t_bc1_us.load() / 1e6);
    s += fmt("barcode 2:             %.1f s\n", g_stats.t_bc2_us.load() / 1e6);
    s += fmt("barcode 3:             %.1f s\n", g_stats.t_bc3_us.load() / 1e6);
    s += fmt("whitelist scoring:     %.1f s\n", g_stats.t_whitelist_us.load() / 1e6);
    s += fmt("Total wall time:       %.2f s\n", wall);
    s += "-----------------------------------------------------\n";
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    std::cout << kBanner << "       preBroCOLI  Version: " << kVersion << "\n" << std::flush;

    Options opt;
    int c, idx = 0;
    while ((c = getopt_long(argc, argv, "q:I:w:x:y:z:o:n:p:f:i:s:c:h",
                            kLongOpts, &idx)) != -1) {
        switch (c) {
            case 'q': opt.mode = optarg; break;
            case 'I': opt.reads = optarg; break;
            case 'w': opt.whitelist = optarg; break;
            case 'x': opt.barcodeX = optarg; break;
            case 'y': opt.barcodeY = optarg; break;
            case 'z': opt.barcodeZ = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 'n': opt.prefix = optarg; break;
            case 'p': opt.threads = std::atoi(optarg); break;
            case 'f': opt.flank_editd = std::atoi(optarg); opt.flank_editd_set = true; break;
            case 'i': opt.trim_barcodes = boolArg(optarg); break;
            case 's': opt.split_by_barcode = boolArg(optarg); break;
            case 'c': opt.mark_chimeric = boolArg(optarg); break;
            case 1001: opt.batch = std::atoi(optarg); break;
            case 1002: opt.keep_unmatched = true; break;
            case 'h': printUsage(argv[0]); return 0;
            default: printUsage(argv[0]); return EXIT_FAILURE;
        }
    }

    // getopt already knows where the options stopped. 1.x counted "params += 2"
    // by hand, which broke as soon as an argument was attached to its flag
    // (-p16) and silently made the tool wait on stdin.
    if (opt.reads.empty()) {
        if (optind < argc) opt.reads = argv[optind];
        else opt.reads = "-";
    }
    if (opt.mode.empty()) {
        std::cerr << "[preBroCOLI] -q/--chemistry is required\n";
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    if (opt.threads < 1) opt.threads = 1;
    if (opt.batch < 1) opt.batch = 1;
    if (!opt.flank_editd_set) opt.flank_editd = (opt.mode == "magicseq") ? 20 : 8;

    LogRedirect log_redirect;

    if (!makeDir(opt.outdir)) return EXIT_FAILURE;

    std::cout << "*****\n"
              << "Chemistry:       " << opt.mode << '\n'
              << "Reads:           " << (opt.reads == "-" ? "<stdin>" : opt.reads) << '\n'
              << "Whitelist:       " << (opt.whitelist.empty() ? "-" : opt.whitelist) << '\n'
              << "Output dir:      " << opt.outdir << '\n'
              << "Prefix:          " << opt.prefix << '\n'
              << "Threads:         " << opt.threads << '\n'
              << "Flank max editd: " << opt.flank_editd << '\n'
              << "Trim barcodes:   " << (opt.trim_barcodes ? "yes" : "no") << '\n'
              << "Keep unmatched:  " << (opt.keep_unmatched ? "yes" : "no") << '\n'
              << "*****\n";

    WhiteListMap whitelist;
    BarcodeSet bcX, bcY, bcZ;
    BarcodeUMIindex buindex;
    loadChemistry(opt, whitelist, bcX, bcY, bcZ, buindex);

    const Pattern pattern = searchPattern(opt.mode);
    for (const auto& p : pattern) std::cout << "  " << p.first << ": " << p.second << "\n";

    if (opt.split_by_barcode && bcX.size() > 50) {
        std::cout << "[preBroCOLI] too many barcodes to split into separate files ("
                  << bcX.size() << " > 50), writing one file instead\n";
        opt.split_by_barcode = false;
    }

    const auto t0 = std::chrono::steady_clock::now();

    std::unique_ptr<SeqReader> reader;
    try {
        reader.reset(new SeqReader(opt.reads));
    } catch (const std::exception& e) {
        std::cerr << "[preBroCOLI] " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    SeqRecord rec;
    // The reader needs one record before isFastq() is meaningful.
    bool have_first = false;
    try {
        have_first = reader->next(rec);
    } catch (const std::exception& e) {
        std::cerr << "[preBroCOLI] " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    if (!have_first) {
        std::cerr << "[preBroCOLI] input has no reads\n";
        return EXIT_FAILURE;
    }
    const bool is_fastq = reader->isFastq();

    OutputWriter out(opt, is_fastq);

    std::cout << "[preBroCOLI] searching for barcodes\n";

    std::vector<std::vector<SearchResult>> batches(static_cast<size_t>(opt.threads));
    std::vector<std::thread> threads(static_cast<size_t>(opt.threads));
    bool eof = false;

    auto take = [&](SearchResult& sr) -> bool {
        if (!have_first) return false;
        sr.read_id = std::move(rec.id);
        sr.line = std::move(rec.seq);
        sr.qual_scores = std::move(rec.qual);
        bump(g_stats.reads_in);
        const long long n = g_stats.reads_in.load();
        if ((n < 1000000 && n % 100000 == 0) || n % 1000000 == 0)
            logLine(fmt("[preBroCOLI]   %.1f million reads read", n / 1e6));
        have_first = reader->next(rec);
        return true;
    };

    try {
        while (!eof) {
            // Every batch has to be cleared before filling, not just the ones
            // this round gets round to. When the input runs out part-way
            // through the loop below it breaks out early, and any batch after
            // that point would still be holding the previous round's reads -
            // which were then searched and written a second time, producing
            // byte-identical duplicate records.
            for (auto& b : batches) b.clear();

            size_t filled = 0;
            for (size_t t = 0; t < batches.size(); ++t) {
                while (static_cast<int>(batches[t].size()) < opt.batch) {
                    SearchResult sr;
                    if (!take(sr)) { eof = true; break; }
                    batches[t].push_back(std::move(sr));
                    ++filled;
                }
                if (eof) break;
            }
            if (filled == 0) break;

            for (size_t t = 0; t < batches.size(); ++t)
                if (!batches[t].empty())
                    threads[t] = std::thread(search_reads, std::cref(opt.mode),
                                             std::ref(batches[t]), std::cref(bcX),
                                             std::cref(bcY), std::cref(bcZ),
                                             opt.flank_editd, std::cref(pattern),
                                             std::cref(whitelist), std::ref(buindex));
            for (size_t t = 0; t < batches.size(); ++t)
                if (!batches[t].empty()) threads[t].join();

            // Writing happens on this thread only, in input order, so the
            // output is reproducible regardless of thread count.
            for (auto& batch : batches)
                for (auto& r : batch) out.write(r);
        }
    } catch (const std::exception& e) {
        std::cerr << "[preBroCOLI] " << e.what() << "\n";
        out.close();
        return EXIT_FAILURE;
    }

    out.close();

    // barcode -> read count, sorted by count desc then barcode asc.
    {
        std::vector<std::pair<std::string, long long>> v(out.barcodeCounts().begin(),
                                                         out.barcodeCounts().end());
        std::sort(v.begin(), v.end(), [](const auto& l, const auto& r) {
            if (l.second != r.second) return l.second > r.second;
            return l.first < r.first;
        });
        LineSink bc(joinPath(opt.outdir, opt.prefix + "_barcode_counts.tsv"));
        bc.append("Barcode\tReads\n");
        for (const auto& p : v) bc.append(p.first + "\t" + std::to_string(p.second) + "\n");
        bc.close();
        std::cout << "[preBroCOLI] barcodes seen: " << v.size() << "\n";
    }

    const double wall = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0).count();
    const std::string summary = summaryText(opt, wall);
    std::cout << summary << std::flush;
    {
        std::ofstream f(joinPath(opt.outdir, opt.prefix + "_summary.txt"), std::ios::trunc);
        f << summary;
    }
    return 0;
}
