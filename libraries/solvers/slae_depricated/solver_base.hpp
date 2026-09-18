#pragma once

#include "MPI_utils.hpp"
#include "OMP_utils.hpp"

#include <Eigen/Sparse>
#include <optional>

namespace nonlocal::slae {

template<class T, class I>
class solver_base {
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& _matrix;
    parallel::MPI_ranges _process_rows;
    size_t _threads_count = parallel::threads_count();

public:
    explicit solver_base(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
        : _matrix{matrix}
        , _process_rows{parallel::rows_distribution(matrix.rows())} {}
    virtual ~solver_base() noexcept = default;

    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix() const noexcept {
        return _matrix;
    }

    const parallel::MPI_ranges& processes_ranges() const noexcept {
        return _process_rows;
    }

    std::ranges::iota_view<size_t, size_t> process_rows(const size_t process = parallel::MPI_rank()) const {
        return _process_rows.get(process);
    }

    size_t threads_count() const noexcept {
        return _threads_count;
    }

    virtual void set_threads_count(const size_t threads_count) {
        _threads_count = threads_count;
    }

    virtual Eigen::Matrix<T, Eigen::Dynamic, 1> solve(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0 = std::nullopt) const = 0;
};

}