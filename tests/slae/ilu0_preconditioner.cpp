#include "create_matrix.hpp"

#include <solvers/slae/ilu0_preconditioner.hpp>

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
// [ 2 10   -1 -0.5    0    0    0    0    0    0]
// [ 1  2   10   -1 -0.5    0    0    0    0    0]
// [ 0  1    2   10   -1 -0.5    0    0    0    0]
// [ 0  0    1    2   10   -1 -0.5    0    0    0]
// [ 0  0    0    1    2   10   -1 -0.5    0    0]
// [ 0  0    0    0    1    2   10   -1 -0.5    0]
// [ 0  0    0    0    0    1    2   10   -1 -0.5]
// [ 0  0    0    0    0    0    1    2   10   -1]
// [ 0  0    0    0    0    0    0    1    2   10]
sparse_matrix<T> scalar_matrix() {
    sparse_matrix<T> matrix{10, 10};
    matrix.portrait.shifts  = {0, 3, 7, 12, 17, 22, 27, 32, 37, 41, 44};
    matrix.portrait.indices = {
        0, 1, 2,
        0, 1, 2, 3,
        0, 1, 2, 3, 4,
        1, 2, 3, 4, 5,
        2, 3, 4, 5, 6,
        3, 4, 5, 6, 7,
        4, 5, 6, 7, 8,
        5, 6, 7, 8, 9,
        6, 7, 8, 9,
        7, 8, 9
    };
    matrix.values = {
                10., -1., -0.5,
            2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,
        1., 2., 10., -1.,
        1., 2., 10.
    };
    return matrix;
}

// 10x10 block matrix with 2x2 blocks
// [10 -1| -0.5    0|    0    0|    0    0|    0    0]
// [ 2 10|   -1 -0.5|    0    0|    0    0|    0    0]
// ---------------------------------------------------
// [ 1  2|   10   -1| -0.5    0|    0    0|    0    0]
// [ 0  1|    2   10|   -1 -0.5|    0    0|    0    0]
// ---------------------------------------------------
// [ 0  0|    1    2|   10   -1| -0.5    0|    0    0]
// [ 0  0|    0    1|    2   10|   -1 -0.5|    0    0]
// ---------------------------------------------------
// [ 0  0|    0    0|    1    2|   10   -1| -0.5    0]
// [ 0  0|    0    0|    0    1|    2   10|   -1 -0.5]
// ---------------------------------------------------
// [ 0  0|    0    0|    0    0|    1    2|   10   -1]
// [ 0  0|    0    0|    0    0|    0    1|    2   10]
sparse_matrix<square_matrix<T, 2>> block_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{5, 5};
    matrix.portrait.shifts  = {0, 2, 5, 8, 11, 13};
    matrix.portrait.indices = {
        0, 1, 
        0, 1, 2, 
        1, 2, 3, 
        2, 3, 4, 
        3, 4
    };
    matrix.values = {
                                             square_matrix<T, 2>{10., -1., 2., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{1., 2., 0., 1.}, square_matrix<T, 2>{10., -1., 2., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{1., 2., 0., 1.}, square_matrix<T, 2>{10., -1., 2., 10.}, square_matrix<T, 2>{-0.5, 0., -1., -0.5},
        square_matrix<T, 2>{1., 2., 0., 1.}, square_matrix<T, 2>{10., -1., 2., 10.},
        square_matrix<T, 2>{1., 2., 0., 1.}, square_matrix<T, 2>{10., -1., 2., 10.}
    };
    return matrix;
}

suite<"ilu0_preconditioner"> _ilu0 = [] {
    "scalar_factorization"_test = [] {
        const auto matrix = scalar_matrix();
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> b = matrix * expected;
        const ilu0_preconditioner preconditioner{scalar_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ilu0 scalar exact solve failed, diff=" << diff;

        print_matrix(preconditioner.matrix());
    };

    "block_factorization"_test = [] {
        const auto matrix = block_matrix();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = matrix * expected;
        const ilu0_preconditioner preconditioner{block_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ilu0 block exact solve failed, diff=" << diff;

        print_matrix(preconditioner.matrix());
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { ilu0_preconditioner{sparse_matrix<T>{3, 4}}; })) <<
            "ilu0 preconditioner must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const ilu0_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const std::vector<T> rhs{1., 2.};
        expect(throws<std::invalid_argument>([&preconditioner, &rhs] { preconditioner.solve(rhs); })) <<
            "ilu0 preconditioner must throw for rhs vector of wrong size";
    };
};

}