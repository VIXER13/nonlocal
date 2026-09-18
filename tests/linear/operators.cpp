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

    "sparse_matrix_scalar_multiplication_block"_test = [] {
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait.shifts = {0, 2, 4, 5};
        block_matrix.portrait.indices = {0, 1, 0, 2, 1};
        block_matrix.values = {{1.0, 2.0, 3.0, 4.0}, 
                               {5.0, 6.0, 7.0, 8.0},
                               {9.0, 10.0, 11.0, 12.0},
                               {13.0, 14.0, 15.0, 16.0},
                               {17.0, 18.0, 19.0, 20.0}};

        block_matrix *= 2.0;
        sparse_matrix<square_matrix<T, 2>> expected_block_matrix{3, 3};
        expected_block_matrix.portrait.shifts = {0, 2, 4, 5};
        expected_block_matrix.portrait.indices = {0, 1, 0, 2, 1};
        expected_block_matrix.values = {{2.0, 4.0, 6.0, 8.0}, 
                                        {10.0, 12.0, 14.0, 16.0},
                                        {18.0, 20.0, 22.0, 24.0},
                                        {26.0, 28.0, 30.0, 32.0},
                                        {34.0, 36.0, 38.0, 40.0}};
        for(const size_t block : std::ranges::iota_view{0u, block_matrix.values.size()})
            for(const size_t row : std::ranges::iota_view{0u, 2u})
                for(const size_t col : std::ranges::iota_view{0u, 2u})
                    expect(approx(block_matrix.values[block][row][col], expected_block_matrix.values[block][row][col], Epsilon)) << 
                        " at block " << block << ", element (" << row << ", " << col << ")";

        block_matrix /= 2.0;
        expected_block_matrix.values = {{1.0, 2.0, 3.0, 4.0}, 
                                        {5.0, 6.0, 7.0, 8.0},
                                        {9.0, 10.0, 11.0, 12.0},
                                        {13.0, 14.0, 15.0, 16.0},
                                        {17.0, 18.0, 19.0, 20.0}};
         for(const size_t block : std::ranges::iota_view{0u, block_matrix.values.size()})
            for(const size_t row : std::ranges::iota_view{0u, 2u})
                for(const size_t col : std::ranges::iota_view{0u, 2u})
                    expect(approx(block_matrix.values[block][row][col], expected_block_matrix.values[block][row][col], Epsilon)) << 
                        " at block " << block << ", element (" << row << ", " << col << ")";
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
        // [1 4|1 2|3 0] * [3] = [35]
        // [3 0|6 7|2 3]   [4]   [77]
        // [-----------]
        // [0 0|4 5|0 0]   [5]   [32]
        // [0 0|6 7|0 0]   [6]   [46]
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait.shifts = {0, 2, 5, 6};
        block_matrix.portrait.indices = {0, 1, 0, 1, 2, 1};
        block_matrix.values = {{1.0, 2.0, 2.0, 2.0}, 
                               {1.0, 4.0, 0.0, 3.0}, 
                               {1.0, 4.0, 3.0, 0.0},
                               {1.0, 2.0, 6.0, 7.0},
                               {3.0, 0.0, 2.0, 3.0},
                               {4.0, 5.0, 6.0, 7.0}};
        const std::vector<std::array<T, 2>> block_vector = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto block_result = block_matrix * block_vector;
        const std::vector<std::array<T, 2>> Expected_Block_Result = {{24.0, 18.0}, {35.0, 77.0}, {32.0, 46.0}};
        expect(eq(block_result.size(), 3));
        for(const size_t row : std::ranges::iota_view{0zu, block_result.size()})
            for(const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(block_result[row][col], Expected_Block_Result[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "sparse_matrix_upper_self_adjoint_vector_multiplication_block"_test = [] {
        // [1 2|1 4|0 0]   [1]   [24]
        // [2 2|0 3|0 0]   [2]   [18]
        // [-----------]
        // [1 0|1 2|3 0] * [3] = [27]
        // [4 3|2 7|2 3]   [4]   [72]
        // [-----------]
        // [0 0|3 2|0 0]   [5]   [17]
        // [0 0|0 3|0 0]   [6]   [12]
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait.shifts = {0, 2, 5, 6};
        block_matrix.portrait.indices = {0, 1, 0, 1, 2, 1};
        block_matrix.values = {{1.0, 2.0, 2.0, 2.0}, 
                               {1.0, 4.0, 0.0, 3.0}, 
                               {1.0, 4.0, 3.0, 0.0},
                               {1.0, 2.0, 6.0, 7.0},
                               {3.0, 0.0, 2.0, 3.0},
                               {4.0, 5.0, 6.0, 7.0}};
        const std::vector<std::array<T, 2>> block_vector = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto result = block_matrix.self_adjoint<matrix_part::Upper>() * block_vector;
        const std::vector<std::array<T, 2>> Expected_Result = {{24.0, 18.0}, {27.0, 72.0}, {17.0, 12.0}};
        expect(eq(result.size(), 3));
        for(const size_t row : std::ranges::iota_view{0zu, result.size()})
            for(const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(result[row][col], Expected_Result[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "sparse_matrix_lower_self_adjoint_vector_multiplication_block"_test = [] {
        // [1 2|1 3|0 0]   [1]   [ 20]
        // [2 2|4 0|0 0]   [2]   [ 18]
        // [-----------]
        // [1 4|1 6|4 6] * [3] = [ 92]
        // [3 0|6 7|5 7]   [4]   [116]
        // [-----------]
        // [0 0|4 5|0 0]   [5]   [ 32]
        // [0 0|6 7|0 0]   [6]   [ 46]
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait.shifts = {0, 2, 5, 6};
        block_matrix.portrait.indices = {0, 1, 0, 1, 2, 1};
        block_matrix.values = {{1.0, 2.0, 2.0, 2.0}, 
                               {1.0, 4.0, 0.0, 3.0}, 
                               {1.0, 4.0, 3.0, 0.0},
                               {1.0, 2.0, 6.0, 7.0},
                               {3.0, 0.0, 2.0, 3.0},
                               {4.0, 5.0, 6.0, 7.0}};
        const std::vector<std::array<T, 2>> block_vector = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto result = block_matrix.self_adjoint<matrix_part::Lower>() * block_vector;
        const std::vector<std::array<T, 2>> Expected_Result = {{20.0, 18.0}, {92.0, 116.0}, {32.0, 46.0}};
        expect(eq(result.size(), 3));
        for(const size_t row : std::ranges::iota_view{0zu, result.size()})
            for(const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(result[row][col], Expected_Result[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "sparse_matrix_vector_multiplication_invalid_size"_test = [] {
        const auto matrix = init_sparse_matrix();
        const std::vector<T> invalid_vector = {1.0, 2.0};
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix * invalid_vector; }));
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix.self_adjoint<matrix_part::Upper>() * invalid_vector; }));
        expect(throws<std::invalid_argument>([&matrix, &invalid_vector] { matrix.self_adjoint<matrix_part::Lower>() * invalid_vector; }));
    };

    "sparse_matrix_addition"_test = [] {
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

    "sparse_matrix_addition_block"_test = [] {
        // [1 2|0 0|0 0]   [0 2|0 0|0 1]   [1 4|0 0|0 1]
        // [2 2|0 0|0 0]   [2 0|0 0|0 2] = [4 2|0 0|0 2]
        // [-----------]   [-----------]   [-----------]
        // [2 0|1 2|0 0] + [2 0|0 4|0 0]   [4 0|1 6|0 0]
        // [0 5|2 3|0 0]   [1 0|4 5|0 0]   [1 5|6 8|0 0]
        // [-----------]   [-----------]   [-----------]
        // [0 0|1 5|0 0]   [0 0|8 2|1 2]   [0 0|9 7|1 2]
        // [0 0|6 7|0 0]   [0 0|1 2|3 4]   [0 0|7 9|3 4]
        sparse_matrix<square_matrix<T, 2>> block_matrix_a{3, 3};
        block_matrix_a.portrait.shifts = {0, 1, 3, 4};
        block_matrix_a.portrait.indices = {0, 0, 1, 1};
        block_matrix_a.values = {{1.0, 2.0, 2.0, 2.0}, 
                                 {2.0, 0.0, 0.0, 5.0}, 
                                 {1.0, 2.0, 2.0, 3.0},
                                 {1.0, 5.0, 6.0, 7.0}};
        sparse_matrix<square_matrix<T, 2>> block_matrix_b{3, 3};
        block_matrix_b.portrait.shifts = {0, 2, 4, 6};
        block_matrix_b.portrait.indices = {0, 2, 0, 1, 1, 2};
        block_matrix_b.values = {{0.0, 2.0, 2.0, 0.0}, 
                                 {0.0, 1.0, 0.0, 2.0}, 
                                 {2.0, 0.0, 1.0, 0.0},
                                 {0.0, 4.0, 4.0, 5.0},
                                 {8.0, 2.0, 1.0, 2.0},
                                 {1.0, 2.0, 3.0, 4.0}};

        block_matrix_a += block_matrix_b;
        expect(eq(block_matrix_a.rows(), 3));
        expect(eq(block_matrix_a.cols(), 3));
        expect(eq(block_matrix_a.non_zeros(), 6));
        expect(approx(block_matrix_a(0, 0)[0][0], 1.0, Epsilon));
        expect(approx(block_matrix_a(0, 0)[0][1], 4.0, Epsilon));
        expect(approx(block_matrix_a(0, 0)[1][0], 4.0, Epsilon));
        expect(approx(block_matrix_a(0, 0)[1][1], 2.0, Epsilon));
        expect(!block_matrix_a.portrait.contains(0, 1));
        expect(approx(block_matrix_a(0, 2)[0][0], 0.0, Epsilon));
        expect(approx(block_matrix_a(0, 2)[0][1], 1.0, Epsilon));
        expect(approx(block_matrix_a(0, 2)[1][0], 0.0, Epsilon));
        expect(approx(block_matrix_a(0, 2)[1][1], 2.0, Epsilon));
        expect(approx(block_matrix_a(1, 0)[0][0], 4.0, Epsilon));
        expect(approx(block_matrix_a(1, 0)[0][1], 0.0, Epsilon));
        expect(approx(block_matrix_a(1, 0)[1][0], 1.0, Epsilon));
        expect(approx(block_matrix_a(1, 0)[1][1], 5.0, Epsilon));
        expect(approx(block_matrix_a(1, 1)[0][0], 1.0, Epsilon));
        expect(approx(block_matrix_a(1, 1)[0][1], 6.0, Epsilon));
        expect(approx(block_matrix_a(1, 1)[1][0], 6.0, Epsilon));
        expect(approx(block_matrix_a(1, 1)[1][1], 8.0, Epsilon));
        expect(!block_matrix_a.portrait.contains(1, 2));
        expect(!block_matrix_a.portrait.contains(2, 0));
        expect(approx(block_matrix_a(2, 1)[0][0], 9.0, Epsilon));
        expect(approx(block_matrix_a(2, 1)[0][1], 7.0, Epsilon));
        expect(approx(block_matrix_a(2, 1)[1][0], 7.0, Epsilon));
        expect(approx(block_matrix_a(2, 1)[1][1], 9.0, Epsilon));
        expect(approx(block_matrix_a(2, 2)[0][0], 1.0, Epsilon));
        expect(approx(block_matrix_a(2, 2)[0][1], 2.0, Epsilon));
        expect(approx(block_matrix_a(2, 2)[1][0], 3.0, Epsilon));
        expect(approx(block_matrix_a(2, 2)[1][1], 4.0, Epsilon));
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