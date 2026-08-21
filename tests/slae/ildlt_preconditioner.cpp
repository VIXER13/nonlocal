#include "create_matrix.hpp"

#include <solvers/slae/ildlt_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;
    
constexpr auto Inf = metamath::constants::Infinity<size_t>;

void print_matrix(const sparse_matrix<T>& matrix) {
    std::cerr << "Matrix: " << std::endl;
    for(size_t row = 0; row < matrix.rows(); ++row) {
        for(const size_t col : matrix.portrait.indices_range(row))
            std::cerr << matrix(row, col) << " ";
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

void print_matrix(const sparse_matrix<square_matrix<T, 2>>& matrix) {
    std::cerr << "Matrix: " << std::endl;
    for(size_t row = 0; row < matrix.rows(); ++row) {
        for(const size_t dof : std::ranges::iota_view{0zu, 2zu}) {
            for(const size_t col : matrix.portrait.indices_range(row)) {
                const auto& val = matrix(row, col);
                std::cerr << val[dof][0] << " " << val[dof][1] << " ";
            }
            std::cerr << std::endl;
        }
    }
    std::cerr << std::endl;
}

// 10x10 general matrix
// [10 -1 -0.5    0    0    0    0    0    0    0]
// [ 0 10   -1 -0.5    0    0    0    0    0    0]
// [ 0  0   10   -1 -0.5    0    0    0    0    0]
// [ 0  0    0   10   -1 -0.5    0    0    0    0]
// [ 0  0    0    0   10   -1 -0.5    0    0    0]
// [ 0  0    0    0    0   10   -1 -0.5    0    0]
// [ 0  0    0    0    0    0   10   -1 -0.5    0]
// [ 0  0    0    0    0    0    0   10   -1 -0.5]
// [ 0  0    0    0    0    0    0    0   10   -1]
// [ 0  0    0    0    0    0    0    0    0   10]
sparse_matrix<T> scalar_matrix() {
    sparse_matrix<T> matrix{10, 10};
    matrix.portrait.shifts  = {0, 3, 6, 9, 12, 15, 18, 21, 24, 26, 27};
    matrix.portrait.indices = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4,
        3, 4, 5,
        4, 5, 6,
        5, 6, 7,
        6, 7, 8,
        7, 8, 9,
        8, 9,
        9
    };
    matrix.values = {
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1., -0.5,
        10., -1.,
        10.
    };
    return matrix;
}

// 10x10 block matrix with 2x2 blocks
// [10 -1| -0.5    0|    0    0|    0    0|    0    0]
// [-1 10|   -1 -0.5|    0    0|    0    0|    0    0]
// ---------------------------------------------------
// [ 0  0|   10   -1| -0.5    0|    0    0|    0    0]
// [ 0  0|   -1   10|   -1 -0.5|    0    0|    0    0]
// ---------------------------------------------------
// [ 0  0|    0    0|   10   -1| -0.5    0|    0    0]
// [ 0  0|    0    0|   -1   10|   -1 -0.5|    0    0]
// ---------------------------------------------------
// [ 0  0|    0    0|    0    0|   10   -1| -0.5    0]
// [ 0  0|    0    0|    0    0|   -1   10|   -1 -0.5]
// ---------------------------------------------------
// [ 0  0|    0    0|    0    0|    0    0|   10   -1]
// [ 0  0|    0    0|    0    0|    0    0|   -1   10]
sparse_matrix<square_matrix<T, 2>> block_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{5, 5};
    matrix.portrait.shifts  = {0, 2, 4, 6, 8, 9};
    matrix.portrait.indices = {
        0, 1, 
        1, 2, 
        2, 3, 
        3, 4, 
        4
    };
    matrix.values = {
        square_matrix<T, 2>{10., -1., -1., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{10., -1., -1., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{10., -1., -1., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{10., -1., -1., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{10., -1., -1., 10.}
    };
    return matrix;
}

suite<"ildlt_preconditioner"> _ildlt = [] {
    "scalar_factorization"_test = [] {
        // ILU0/ILDLT is exact for banded matrices: (L D L^T)^{-1} A x = x.
        const auto matrix = scalar_matrix();
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{scalar_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 8.9e-16;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt scalar exact solve failed, diff=" << diff;

        print_matrix(preconditioner.matrix());
    };

    "block_factorization"_test = [] {
        const auto matrix = block_matrix();
        validate_sparse_matrix(matrix);

        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{block_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt block exact solve failed, diff=" << diff;

        print_matrix(preconditioner.matrix());
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { ildlt_preconditioner{sparse_matrix<T>{3, 4}}; })) <<
            "ildlt preconditioner must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const ildlt_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const std::vector<T> rhs{1., 2.};
        expect(throws<std::invalid_argument>([&preconditioner, &rhs] { preconditioner.solve(rhs); })) <<
            "ildlt preconditioner must throw for rhs vector of wrong size";
    };
};

}