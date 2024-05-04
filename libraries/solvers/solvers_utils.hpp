#pragma once

#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <ranges>
#include <vector>
#include <variant>

namespace nonlocal::utils {

class nodes_sequence : public std::variant<std::ranges::iota_view<size_t, size_t>,
                                           std::vector<size_t>> {
    using _base = std::variant<
        std::ranges::iota_view<size_t, size_t>,
        std::vector<size_t>
    >;

public:
    using _base::_base;
    using _base::operator=;

    constexpr bool is_range() const noexcept {
        return std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(*this);
    }

    constexpr bool is_vector() const noexcept {
        return std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(*this);
    }

    constexpr std::ranges::iota_view<size_t, size_t> as_range() const {
        return std::get<std::ranges::iota_view<size_t, size_t>>(*this);
    }

    constexpr std::vector<size_t>& as_vector() {
        return std::get<std::vector<size_t>>(*this);
    }

    constexpr const std::vector<size_t>& as_vector() const {
        return std::get<std::vector<size_t>>(*this);
    }
};

// using nodes_sequence = std::variant<
//     std::ranges::iota_view<size_t, size_t>,
//     std::vector<size_t>
// >;

template<class T, class I>
void accumulate_shifts(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
    for(const size_t i : std::ranges::iota_view{0u, size_t(K.rows())})
        K.outerIndexPtr()[i + 1] += K.outerIndexPtr()[i];
}

template<class T, class I>
void clear_matrix_coefficients(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
    for(const size_t i : std::ranges::iota_view{size_t{0}, size_t(K.nonZeros())})
        K.valuePtr()[i] = T{0};
}

template<class Sequence, class T, class I>
void clear_matrix_rows(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, 
                       const Sequence& rows) {
    for(const size_t row : rows)
        for(const size_t shift : std::ranges::iota_view{matrix.outerIndexPtr()[row], matrix.outerIndexPtr()[row + 1]})
            matrix.valuePtr()[shift] = T{0};
}

template<class T, class I>
void clear_matrix_rows(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, 
                       const nodes_sequence& rows) {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(rows))
        clear_matrix_rows(matrix, std::get<std::ranges::iota_view<size_t, size_t>>(rows));
    else
        clear_matrix_rows(matrix, std::get<std::vector<size_t>>(rows));
}

template<class T, class I>
void multiply_matrix_rows(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                          const std::ranges::iota_view<size_t, size_t> rows,
                          const T value) {
    for(const size_t row : rows)
        for(const size_t shift : std::ranges::iota_view{matrix.outerIndexPtr()[row], matrix.outerIndexPtr()[row + 1]})
            matrix().valuePtr()[shift] *= value;
}

template<class T, class I>
void allocate_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
    K.data().resize(K.outerIndexPtr()[K.rows()]);
    for(const size_t i : std::ranges::iota_view{size_t{0}, size_t(K.nonZeros())}) {
        K.innerIndexPtr()[i] = 0;
        K.valuePtr()[i] = T{0};
    }
}

template<class T, class I>
void sort_indices(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
#pragma omp parallel for default(none) shared(K) schedule(dynamic)
    for(size_t i = 0; i < K.rows(); ++i)
        std::sort(&K.innerIndexPtr()[K.outerIndexPtr()[i]], &K.innerIndexPtr()[K.outerIndexPtr()[i+1]]);
}

}