#pragma once

#include "iterative_solver_base.hpp"

#include <parallel/OMP_utils.hpp>
#include <parallel/uniform_ranges.hpp>
#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/utils/operators.hpp>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace nonlocal::slae {

template<class T, std::integral I = uint32_t, std::integral J = size_t>
class matrix_vector_production : public iterative_solver_base<T, I, J> {
    using _base = iterative_solver_base<T, I, J>;

    parallel::OMP_ranges _thread_rows;

public:
    using typename _base::entity_t;
    using _base::threads_count;
    using _base::matrix;

protected:
    void production(std::vector<entity_t>& result, const std::vector<entity_t>& vector) const {
    #pragma omp parallel num_threads(threads_count())
    {
    #ifdef _OPENMP
        const int thread = omp_get_thread_num();
    #else
        const int thread = 0;
    #endif
        for(const size_t row : _thread_rows.get(thread)) {
            result[row] = {};
            using namespace metamath::linear;
            using namespace metamath::operators;
            for(const size_t shift : matrix().portrait.shifts_range(row))
                result[row] += matrix().values[shift] * vector[matrix().portrait.indices[shift]];
        }
    }
    }

public:
    explicit matrix_vector_production(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _base{matrix}
        , _thread_rows{parallel::uniform_ranges(matrix.portrait.shifts, threads_count())} {}
};

template<class T, std::integral I = uint32_t, std::integral J = size_t>
class symmetric_matrix_vector_product : public iterative_solver_base<T, I, J> {
    using _base = iterative_solver_base<T, I, J>;

public:
    using typename _base::entity_t;
    using _base::threads_count;
    using _base::matrix;

private:
    parallel::OMP_ranges _thread_rows;
    mutable metamath::linear::matrix<entity_t> _threaded_product;

protected:
    void production(std::vector<entity_t>& result, const std::vector<entity_t>& vector) const {
    #pragma omp parallel num_threads(threads_count())
    {
    #ifdef _OPENMP
        const int thread = omp_get_thread_num();
    #else
        const int thread = 0;
    #endif
        for(const size_t row : _thread_rows.get(thread))
            for(const size_t shift : matrix().portrait.shifts_range(row)) {
                using namespace metamath::linear;
                using namespace metamath::operators;
                const size_t col = matrix().portrait.indices[shift];
                _threaded_product(thread, row) += matrix().values[shift] * vector[col];
                if (col != row)
                    _threaded_product(thread, col) += matrix().values[shift] * vector[row];
            }
    }
        for(const size_t row : std::ranges::iota_view{0u, result.size()}) {
            result[row] = {};
            for(const size_t thread : std::ranges::iota_view{0zu, threads_count()})
                result[row] += _threaded_product(thread, row);
        }
    }

public:
    explicit symmetric_matrix_vector_product(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _base{matrix}
        , _thread_rows{parallel::uniform_ranges(matrix.portrait.shifts, threads_count())}
        , _threaded_product{threads_count(), matrix.cols()} {}
};

}