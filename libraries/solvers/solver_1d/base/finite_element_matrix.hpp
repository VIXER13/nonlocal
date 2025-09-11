#pragma once

#include <solvers/base/utils.hpp>

#include <Eigen/Sparse>

#include <array>
#include <unordered_map>

namespace nonlocal::solver_1d {

template<class T, class I>
struct finite_element_matrix_1d final {
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> inner;
    std::array<std::unordered_map<size_t, T>, 2> bound;

    void clear() {
        inner = {};
        bound = {};
    }

    void set_zero(const std::optional<utils::nodes_sequence>& rows_sequence = std::nullopt) {
        if (!rows_sequence) {
            bound.front().clear();
            bound.back().clear();
            std::memset(inner.valuePtr(), '\0', sizeof(T) * size_t(inner.nonZeros()));
        } 
        else {
            utils::iterate(*rows_sequence, [this](const size_t row) {
                const size_t bytes = sizeof(T) * (inner.outerIndexPtr()[row + 1] - inner.outerIndexPtr()[row]);
                std::memset(inner.valuePtr() + inner.outerIndexPtr()[row], '\0', bytes);
                for (auto& boundary : bound)
                    if (const auto it = boundary.find(row); it != boundary.cend())
                        boundary.erase(it);
            });
        }
    }
};

}