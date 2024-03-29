#include "OMP_utils.hpp"

#include "init_uniform_ranges.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace parallel_utils {

// GCC implementation have bug: omp_get_num_threads() return 1 in nonparallel sections
int threads_count() {
    int threads = 1;
#ifdef _OPENMP
    #pragma omp parallel default(none) shared(threads)
    {
        threads = omp_get_num_threads();
    }
#endif
    return threads > 1 ? threads : 1;
}

OMP_ranges::OMP_ranges(const size_t size, const size_t threads)
    : _ranges{init_uniform_ranges(size, threads)} {}

std::ranges::iota_view<size_t, size_t> OMP_ranges::get(const size_t thread) const {
    return _ranges[thread];
}

void OMP_ranges::set(const std::ranges::iota_view<size_t, size_t> range, const size_t thread) {
    _ranges[thread] = range;
}

}