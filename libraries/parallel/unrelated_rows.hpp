#ifndef PARALLEL_UTILS_UNRELATED_ROWS_HPP
#define PARALLEL_UTILS_UNRELATED_ROWS_HPP

#include <Eigen/Sparse>

#include <cstddef>
#include <ranges>
#include <vector>

namespace parallel {

template<class I>
class unrelated_rows final {
    template<class T>
    static bool has_related_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                 const std::vector<bool>& is_related, 
                                 const std::ranges::iota_view<I, I> indices);

public:
    std::vector<I> shifts; // Independent rows number in a rows sequence below.
    std::vector<I> rows; // The rows sequence in which they must be calculated so that they are independent in parallel code.
                         // The rows count is the same as the rows number in the matrix.

    constexpr unrelated_rows() noexcept = default;
    template<class T>
    explicit unrelated_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);

    template<class T>
    void init(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);

    void clear() noexcept;
};

template<class I>
template<class T>
unrelated_rows<I>::unrelated_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    init(matrix);
}

template<class I>
template<class T>
bool unrelated_rows<I>::has_related_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                         const std::vector<bool>& is_related, 
                                         const std::ranges::iota_view<I, I> indices) {
    for(const I i : indices)
        if (is_related[matrix.innerIndexPtr()[i]])
            return true;
    return false;
}

template<class I>
template<class T>
void unrelated_rows<I>::init(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    clear();
    shifts.push_back(0);
    rows.reserve(matrix.rows());
    std::vector<bool> is_added(matrix.cols(), false);   // .cols() used firstly, because we believe that the matrix is square,
    std::vector<bool> is_related(matrix.cols(), false); // and secondly, because of MPI, where the matrix is cut row by row
    for(const I start_row : std::ranges::iota_view{I{0}, matrix.rows()}) {
        if (is_added[start_row])
            continue;
        I& shift = shifts.emplace_back(shifts.back());
        for(const I row : std::ranges::iota_view{start_row, matrix.rows()}) {
            if (is_added[row])
                continue;
            if (const auto inner_indices = std::ranges::iota_view{matrix.outerIndexPtr()[row], matrix.outerIndexPtr()[row + 1]};
                !has_related_rows(matrix, is_related, inner_indices)) {
                for(const size_t i : inner_indices)
                    is_related[matrix.innerIndexPtr()[i]] = true;
                is_added[row] = true;
                rows.push_back(row);
                ++shift;
            }
        }
        std::fill(is_related.begin(), is_related.end(), false);
    }
}

template<class I>
void unrelated_rows<I>::clear() noexcept {
    shifts.clear();
    rows.clear();
}

}

#endif