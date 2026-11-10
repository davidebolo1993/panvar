#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace panvar {

// Run fn(i) for i in [0, n) across worker threads (work-stealing on an atomic
// counter). threads==0 -> hardware concurrency. Results must be written to
// distinct indices per i (no shared mutation) for thread safety + determinism.
// An exception thrown by fn stops the hand-out of further indices and is rethrown
// here once every worker has finished, as it would be with one thread: escaping a
// worker it would call std::terminate, past every caller's cleanup. When several
// indices fail before the hand-out stops, the lowest of them is rethrown.
inline void run_parallel(std::size_t n, std::size_t threads,
                         const std::function<void(std::size_t)>& fn) {
    std::size_t nthreads = threads != 0 ? threads : std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    if (n == 0) return;
    nthreads = std::min(nthreads, n);
    if (nthreads <= 1) {
        for (std::size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr first_error;
    std::size_t first_error_index = n;
    std::mutex error_mutex;
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    // A thread that cannot be started throws from emplace_back; the ones already running must be
    // joined before that propagates, since destroying a joinable std::thread calls std::terminate.
    try {
        for (std::size_t t = 0; t < nthreads; ++t) {
            pool.emplace_back([&] {
                for (;;) {
                    if (failed.load()) break;
                    const std::size_t i = next.fetch_add(1);
                    if (i >= n) break;
                    try {
                        fn(i);
                    } catch (...) {
                        std::lock_guard<std::mutex> lock(error_mutex);
                        if (!first_error || i < first_error_index) {
                            first_error = std::current_exception();
                            first_error_index = i;
                        }
                        failed.store(true);
                        break;
                    }
                }
            });
        }
    } catch (...) {
        failed.store(true);
        for (std::thread& th : pool) th.join();
        throw;
    }
    for (std::thread& th : pool) th.join();
    if (first_error) std::rethrow_exception(first_error);
}

} // namespace panvar
