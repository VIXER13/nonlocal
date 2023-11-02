#ifndef NONLOCAL_SLAE_SOLVER_BASE_HPP
#define NONLOCAL_SLAE_SOLVER_BASE_HPP

#include "OMP_utils.hpp"

#include <Eigen/Sparse>
#include <optional>

namespace nonlocal::slae {

template<class T, class I>
class solver_base {
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& _matrix;
    size_t _threads_count = parallel_utils::threads_count();

public:
    explicit solver_base(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
        : _matrix{matrix} {}
    virtual ~solver_base() noexcept = default;

    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix() const noexcept {
        return _matrix;
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

#endif