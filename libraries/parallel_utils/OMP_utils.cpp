#include "OMP_utils.hpp"

#include <omp.h>

namespace parallel_utils {

// GCC implementation have bug: omp_get_num_threads() return 1 in nonparallel sections
int threads_count() {
    int threads = 1;
#pragma omp parallel default(none) shared(threads)
{
    threads = omp_get_num_threads();
}
    return threads > 1 ? threads : 1;
}

}