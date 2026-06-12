#include <metamath/linear/linear.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;

using T = double;

constexpr T Epsilon = std::numeric_limits<T>::epsilon();

sparse_matrix<T> init_sparse_matrix() {
    // [1 0 3]
    // [0 2 4]
    // [5 0 0]
    sparse_matrix<T> matrix{3, 3};
    matrix.portrait.shifts = {0, 2, 4, 5};
    matrix.portrait.indices = {0, 2, 1, 2, 0};
    matrix.values = {1.0, 3.0, 2.0, 4.0, 5.0};
    return matrix;
}

std::pair<sparse_matrix<T>, std::vector<T>> init_sparse_matrix_and_vector() {
    return {init_sparse_matrix(), {1.0, 2.0, 3.0}};
}

suite<"sparse_matrix_operators"> _ = [] {
    "sparse_matrix_scalar_multiplication"_test = [] {
        sparse_matrix<T> matrix{3, 3};
        matrix.portrait.shifts = {0, 1, 2, 3};
        matrix.portrait.indices = {0, 1, 2};
        matrix.values = {1.0, 2.0, 3.0};

        matrix *= 2.0;
        expect(approx(matrix.values[0], 2.0, Epsilon));
        expect(approx(matrix.values[1], 4.0, Epsilon));
        expect(approx(matrix.values[2], 6.0, Epsilon));

        matrix /= 2.0;
        expect(approx(matrix.values[0], 1.0, Epsilon));
        expect(approx(matrix.values[1], 2.0, Epsilon));
        expect(approx(matrix.values[2], 3.0, Epsilon));
    };

    "sparse_matrix_vector_multiplication"_test = [] {
        // [1 0 3]   [1]   [10]
        // [0 2 4] * [2] = [16]
        // [5 0 0]   [3]   [ 5]
        const auto [matrix, vector] = init_sparse_matrix_and_vector();
        const auto result = matrix * vector;
        expect(eq(result.size(), 3));
        expect(approx(result[0], 10.0, Epsilon));
        expect(approx(result[1], 16.0, Epsilon));
        expect(approx(result[2], 5.0, Epsilon));
    };

    "sparse_matrix_upper_self_adjoint_vector_multiplication"_test = [] {
        // [1 0 3]   [1]   [10]
        // [0 2 4] * [2] = [16]
        // [3 4 0]   [3]   [11]
        const auto [matrix, vector] = init_sparse_matrix_and_vector();
        const auto upper_triangular_result = matrix.self_adjoint<matrix_part::Upper>() * vector;
        expect(eq(upper_triangular_result.size(), 3));
        expect(approx(upper_triangular_result[0], 10.0, Epsilon));
        expect(approx(upper_triangular_result[1], 16.0, Epsilon));
        expect(approx(upper_triangular_result[2], 11.0, Epsilon));
    };

    "sparse_matrix_lower_self_adjoint_vector_multiplication"_test = [] {
        // [1 0 5]   [1]   [16]
        // [0 2 0] * [2] = [ 4]
        // [5 0 0]   [3]   [ 5]
        const auto [matrix, vector] = init_sparse_matrix_and_vector();
        const auto lower_triangular_result = matrix.self_adjoint<matrix_part::Lower>() * vector;
        expect(eq(lower_triangular_result.size(), 3));
        expect(approx(lower_triangular_result[0], 16.0, Epsilon));
        expect(approx(lower_triangular_result[1], 4.0, Epsilon));
        expect(approx(lower_triangular_result[2], 5.0, Epsilon));
    };

    "sparse_matrix_vector_multiplication_block"_test = [] {
        // [1 2|1 4|0 0]   [1]   [24]
        // [2 2|0 3|0 0]   [2]   [18]
        // [-----------]
        // [2 0|0 0|3 0] * [3] = [17]
        // [0 5|0 0|2 3]   [4]   [38]
        // [-----------]
        // [0 0|4 5|0 0]   [5]   [32]
        // [0 0|6 7|0 0]   [6]   [46]
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait.shifts = {0, 2, 4, 5};
        block_matrix.portrait.indices = {0, 1, 0, 2, 1};
        block_matrix.values = {{1.0, 2.0, 2.0, 2.0}, 
                               {1.0, 4.0, 0.0, 3.0}, 
                               {2.0, 0.0, 0.0, 5.0},
                               {3.0, 0.0, 2.0, 3.0},
                               {4.0, 5.0, 6.0, 7.0}};
        const std::vector<std::array<T, 2>> block_vector = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto block_result = block_matrix * block_vector;
        expect(eq(block_result.size(), 3));
        expect(approx(block_result[0][0], 24.0, Epsilon));
        expect(approx(block_result[0][1], 18.0, Epsilon));
        expect(approx(block_result[1][0], 17.0, Epsilon));
        expect(approx(block_result[1][1], 38.0, Epsilon));
        expect(approx(block_result[2][0], 32.0, Epsilon));
        expect(approx(block_result[2][1], 46.0, Epsilon));
    };

    "sparse_matrix_vector_multiplication_invalid_size"_test = [] {
        const auto matrix = init_sparse_matrix();
        const std::vector<T> invalid_vector = {1.0, 2.0};
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix * invalid_vector; }));
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix.self_adjoint<matrix_part::Upper>() * invalid_vector; }));
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix.self_adjoint<matrix_part::Lower>() * invalid_vector; }));
    };

    "sparse_matrix_addition_4x4"_test = [] {
        // [1 0 3 0]   [0 2 0 1]   [1 2 3 1]
        // [0 2 4 0] + [2 0 3 0] = [2 2 7 0]
        // [5 0 0 0]   [0 5 0 2]   [5 5 0 2]
        // [0 0 0 6]   [1 0 4 0]   [1 0 4 6]
        sparse_matrix<T> matrix_a{4, 4};
        matrix_a.portrait.shifts = {0, 2, 4, 5, 6};
        matrix_a.portrait.indices = {0, 2, 1, 2, 0, 3};
        matrix_a.values = {1.0, 3.0, 2.0, 4.0, 5.0, 6.0};

        sparse_matrix<T> matrix_b{4, 4};
        matrix_b.portrait.shifts = {0, 2, 4, 6, 8};
        matrix_b.portrait.indices = {1, 3, 0, 2, 1, 3, 0, 2};
        matrix_b.values = {2.0, 1.0, 2.0, 3.0, 5.0, 2.0, 1.0, 4.0};

        matrix_a += matrix_b;
        expect(eq(matrix_a.rows(), 4));
        expect(eq(matrix_a.cols(), 4));
        expect(eq(matrix_a.non_zeros(), 13));
        expect(approx(matrix_a(0, 0), 1.0, Epsilon));
        expect(approx(matrix_a(0, 1), 2.0, Epsilon));
        expect(approx(matrix_a(0, 2), 3.0, Epsilon));
        expect(approx(matrix_a(0, 3), 1.0, Epsilon));
        expect(approx(matrix_a(1, 0), 2.0, Epsilon));
        expect(approx(matrix_a(1, 1), 2.0, Epsilon));
        expect(approx(matrix_a(1, 2), 7.0, Epsilon));
        expect(!matrix_a.portrait.contains(1, 3));
        expect(approx(matrix_a(2, 0), 5.0, Epsilon));
        expect(approx(matrix_a(2, 1), 5.0, Epsilon));
        expect(!matrix_a.portrait.contains(2, 2));
        expect(approx(matrix_a(2, 3), 2.0, Epsilon));
        expect(approx(matrix_a(3, 0), 1.0, Epsilon));
        expect(!matrix_a.portrait.contains(3, 1));
        expect(approx(matrix_a(3, 2), 4.0, Epsilon));
        expect(approx(matrix_a(3, 3), 6.0, Epsilon));
    };

    "sparse_matrix_addition_invalid_size"_test = [] {
        sparse_matrix<T> matrix_a{3, 3};
        matrix_a.portrait.shifts = {0, 2, 4, 5};
        matrix_a.portrait.indices = {0, 2, 1, 2, 0};
        matrix_a.values = {1.0, 3.0, 2.0, 4.0, 5.0};

        sparse_matrix<T> matrix_b{4, 4};
        matrix_b.portrait.shifts = {0, 2, 4, 6, 8};
        matrix_b.portrait.indices = {1, 3, 0, 2, 1, 3, 0, 2};
        matrix_b.values = {2.0, 1.0, 2.0, 3.0, 5.0, 2.0, 1.0, 4.0};

        expect(throws<std::invalid_argument>([&matrix_a, &matrix_b] { matrix_a += matrix_b; }));
    };
};

}