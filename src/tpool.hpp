// Thin C++ wrapper around htslib's official thread pool (htslib/thread_pool.h).
//
// Using hts_tpool rather than a hand-rolled pool means the same worker threads
// can, if we want, also serve BGZF (de)compression. We deliberately do NOT
// attach the pool to the per-worker BAM readers: stage 1 already runs one
// reader per thread, so enabling htslib's own decompression threads on top
// would both oversubscribe the CPU and risk a dispatch deadlock (a worker
// blocking on a decompression task that needs a worker).
#pragma once

#include <atomic>
#include <exception>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <string>

#include <htslib/thread_pool.h>

namespace brocoli {

class ThreadPool {
public:
    explicit ThreadPool(int nthreads) : n_(nthreads < 1 ? 1 : nthreads) {
        pool_ = hts_tpool_init(n_);
        if (!pool_) throw std::runtime_error("hts_tpool_init failed");
    }
    ~ThreadPool() { if (pool_) hts_tpool_destroy(pool_); }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    hts_tpool* raw() const { return pool_; }
    int size() const { return n_; }

private:
    hts_tpool* pool_ = nullptr;
    int n_;
};

// A queue of fire-and-forget jobs bound to a pool. `submit` blocks once the
// queue is full, which gives us natural backpressure and bounds memory when a
// producer is much faster than the consumers.
class TaskQueue {
public:
    TaskQueue(ThreadPool& pool, int qsize = 0)
        : pool_(pool.raw()) {
        if (qsize <= 0) qsize = pool.size() * 4;
        // in_only = 1: results are discarded, nothing accumulates in the queue.
        q_ = hts_tpool_process_init(pool_, qsize, 1);
        if (!q_) throw std::runtime_error("hts_tpool_process_init failed");
    }
    ~TaskQueue() { if (q_) hts_tpool_process_destroy(q_); }

    TaskQueue(const TaskQueue&) = delete;
    TaskQueue& operator=(const TaskQueue&) = delete;

    void submit(std::function<void()> fn) {
        auto* job = new Job{std::move(fn), this};
        if (hts_tpool_dispatch(pool_, q_, &TaskQueue::trampoline, job) != 0) {
            delete job;
            throw std::runtime_error("hts_tpool_dispatch failed");
        }
    }

    // Wait for every dispatched job to finish, then rethrow the first error.
    void wait() {
        hts_tpool_process_flush(q_);
        std::lock_guard<std::mutex> lk(err_mu_);
        if (!err_.empty()) {
            const std::string e = err_;
            err_.clear();
            throw std::runtime_error(e);
        }
    }

private:
    struct Job {
        std::function<void()> fn;
        TaskQueue* self;
    };

    static void* trampoline(void* arg) {
        auto* job = static_cast<Job*>(arg);
        try {
            job->fn();
        } catch (const std::exception& e) {
            job->self->recordError(e.what());
        } catch (...) {
            job->self->recordError("unknown exception in worker");
        }
        delete job;
        return nullptr;
    }

    void recordError(const std::string& msg) {
        std::lock_guard<std::mutex> lk(err_mu_);
        if (err_.empty()) err_ = msg;
    }

    hts_tpool* pool_ = nullptr;
    hts_tpool_process* q_ = nullptr;
    std::mutex err_mu_;
    std::string err_;
};

}  // namespace brocoli
