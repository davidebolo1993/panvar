#pragma once

#include <atomic>
#include <cstddef>
#include <functional>
#include <thread>
#include <vector>

namespace panvar {

// Run fn(i) for i in [0, n) across worker threads (work-stealing on an atomic
// counter). threads==0 -> hardware concurrency. Results must be written to
// distinct indices per i (no shared mutation) for thread safety + determinism.
inline void run_parallel(std::size_t n, std::size_t threads,
                         const std::function<void(std::size_t)>& fn) {
    std::size_t nthreads = threads != 0 ? threads : std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    if (n > 0) nthreads = std::min(nthreads, n);
    if (nthreads <= 1) {
        for (std::size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (std::size_t t = 0; t < nthreads; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                const std::size_t i = next.fetch_add(1);
                if (i >= n) break;
                fn(i);
            }
        });
    }
    for (std::thread& th : pool) th.join();
}

} // namespace panvar
