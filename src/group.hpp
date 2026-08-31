// Stage 2: merge every sample's record stream and cut the genome into read
// clusters ("groups").
//
// BroCOLI 1.x formed clusters inside each byte chunk, merged the chunks, then
// tried to reconcile clusters across samples by mapping intervals onto
// intervals (Merge_Read_Interval + get_pointers). Because every stage-1 stream
// is now guaranteed to be in genome order, a plain k-way merge gives the same
// answer in one pass and cannot disagree with itself at a boundary.
#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.hpp"
#include "record.hpp"

namespace brocoli {

// One sample's chunk files presented as a single ordered stream.
class SampleStream {
public:
    SampleStream(std::vector<std::string> files) : files_(std::move(files)) { advance(); }

    bool valid() const { return valid_; }
    const ReadRec& rec() const { return cur_; }
    const std::string& blob() const { return blob_; }

    void advance() {
        valid_ = false;
        while (fi_ <= files_.size()) {
            if (!rd_) {
                if (fi_ >= files_.size()) return;
                rd_.reset(new RecordReader(files_[fi_]));
            }
            if (rd_->next(cur_, &blob_)) { valid_ = true; return; }
            rd_.reset();
            ++fi_;
        }
    }

private:
    std::vector<std::string> files_;
    size_t fi_ = 0;
    std::unique_ptr<RecordReader> rd_;
    ReadRec cur_;
    std::string blob_;
    bool valid_ = false;
};

struct MergeResult {
    std::string record_path;
    std::vector<GroupIndexEntry> groups;
    long long n_records = 0;
};

inline MergeResult mergeAndGroup(const std::vector<std::vector<std::string>>& per_sample,
                                 const std::string& out_path) {
    MergeResult res;
    res.record_path = out_path;

    std::vector<std::unique_ptr<SampleStream>> streams;
    streams.reserve(per_sample.size());
    for (const auto& f : per_sample)
        streams.emplace_back(new SampleStream(f));

    RecordWriter writer(out_path);

    GroupIndexEntry cur{};
    bool open = false;
    int gid = 0;

    auto closeGroup = [&] {
        if (!open) return;
        res.groups.push_back(cur);
        open = false;
    };

    while (true) {
        // pick the smallest (tid, beg) among the sample heads
        int best = -1;
        for (size_t i = 0; i < streams.size(); ++i) {
            if (!streams[i]->valid()) continue;
            if (best < 0) { best = static_cast<int>(i); continue; }
            const ReadRec& a = streams[i]->rec();
            const ReadRec& b = streams[best]->rec();
            if (a.tid < b.tid || (a.tid == b.tid && a.beg < b.beg))
                best = static_cast<int>(i);
        }
        if (best < 0) break;

        const ReadRec& r = streams[best]->rec();

        if (!open || r.tid != cur.tid || r.beg > cur.high) {
            closeGroup();
            cur = GroupIndexEntry{};
            cur.offset  = writer.tell();
            cur.tid     = r.tid;
            cur.low     = r.beg;
            cur.high    = r.end;
            cur.n_reads = 0;
            cur.gid     = ++gid;
            open = true;
        }
        cur.high = std::max(cur.high, r.end);
        ++cur.n_reads;
        ++res.n_records;

        writer.appendRaw(streams[best]->blob().data(), streams[best]->blob().size());
        streams[best]->advance();
    }
    closeGroup();
    writer.close();

    std::cout << "[BroCOLI] " << res.n_records << " reads in "
              << res.groups.size() << " clusters\n";
    return res;
}

// ---------------------------------------------------------------- group ----
// Everything the identification/quantification stage needs about one cluster.
// Reads are keyed by name in bulk mode and by "barcode-UMI" in single-cell
// mode, which reproduces the exact-duplicate collapsing of BroCOLI 1.x.
struct GroupData {
    std::string chrom;
    Interval span{{0, 0}};
    std::string gid_str;

    std::unordered_map<std::string, int>          read_file;
    std::unordered_map<std::string, std::string>  read_bc;
    std::unordered_map<std::string, std::string>  read_umi;
    std::unordered_map<std::string, IntervalVec>  read_sjs;
    std::unordered_map<std::string, std::vector<int>> read_sigs;
    std::unordered_map<std::string, Interval>     read_cov;
    std::unordered_map<std::string, Interval>     single_exon;

    size_t size() const { return read_file.size(); }
};

inline GroupData loadGroup(RecordReader& reader, const GroupIndexEntry& e,
                           const std::string& chrom, Mode mode) {
    GroupData g;
    g.chrom   = chrom;
    g.span    = Interval{e.low, e.high};
    g.gid_str = std::to_string(e.gid);

    const size_t n = static_cast<size_t>(e.n_reads);
    g.read_file.reserve(n * 2);
    g.read_cov.reserve(n * 2);

    reader.seek(e.offset);

    std::unordered_map<std::string, int> best_len;
    best_len.reserve(n * 2);

    ReadRec r;
    for (int i = 0; i < e.n_reads; ++i) {
        if (!reader.next(r)) break;

        const std::string key =
            (mode == Mode::SC) ? (r.bc + "-" + r.umi) : r.name;

        auto bl = best_len.find(key);
        if (bl != best_len.end()) {
            ++g_stats.umi_collapsed;
            if (r.len <= bl->second) continue;   // keep the longest representative
            bl->second = r.len;
            g.read_sjs.erase(key);
            g.read_sigs.erase(key);
            g.single_exon.erase(key);
        } else {
            best_len.emplace(key, r.len);
        }

        g.read_file[key] = r.file;
        g.read_cov[key]  = Interval{r.beg, r.end};
        if (mode == Mode::SC) {
            g.read_bc[key]  = r.bc;
            g.read_umi[key] = r.umi;
        }
        if (r.sjs.empty()) {
            g.single_exon[key] = Interval{r.beg, r.end};
        } else {
            g.read_sjs[key]  = r.sjs;
            g.read_sigs[key] = r.sigs;
        }
    }
    g_stats.single_exon += static_cast<long long>(g.single_exon.size());
    return g;
}

}  // namespace brocoli
