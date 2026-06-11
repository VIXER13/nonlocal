#include <metamath/linear/linear.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;

using T = double;
using I = int32_t;

suite<"sparse_matrix"> _ = [] {
    "common"_test = [] {
        // Matrix representation in that test is the following:
        // [1 0 3 0]
        // [0 2 0 4]
        // [5 0 0 0]
        sparse_matrix<T> matrix;

        // Initially matrix is empty
        expect(eq(matrix.rows(), 0));
        expect(eq(matrix.cols(), 0));
        expect(eq(matrix.non_zeros(), 0));

        // Set matrix size
        matrix.portrait().set_size(3, 4);
        expect(eq(matrix.rows(), 3));
        expect(eq(matrix.cols(), 4));

        // Fill shifts vector with non-accumulated values
        matrix.portrait().shifts() = {0, 2, 2, 1};
        // Check that validation fails due to invalid shifts vector
        expect(throws<std::logic_error>([&matrix] { validate_shifts(matrix.portrait().shifts()); }));
        // Accumulate shifts and check that validation passes
        matrix.portrait().accumulate_shifts();
        expect(eq(matrix.non_zeros(), 5));
        expect(nothrow([&matrix] { validate_shifts(matrix.portrait().shifts()); }));

        // Allocate indices and fill them with unsorted column indices
        matrix.portrait().allocate_indices();
        matrix.portrait().indices()[0] = 2;
        matrix.portrait().indices()[1] = 0;
        matrix.portrait().indices()[2] = 3;
        matrix.portrait().indices()[3] = 1;
        matrix.portrait().indices()[4] = 0;
        // Check that validation fails due to unsorted column indices within rows
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));
        matrix.portrait().sort_indices();
        // Check that after sorting validation passes
        expect(nothrow([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Check that validation fails due to values vector size not matching the number of non-zero elements in the portrait
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix(matrix); }));
        matrix.allocate_values();
        // Allocate values and check that validation passes
        expect(nothrow([&matrix] { validate_sparse_matrix(matrix); }));

        // Check that contains() method work correctly
        expect( matrix.portrait().contains(0, 0));
        expect(!matrix.portrait().contains(0, 1));
        expect( matrix.portrait().contains(0, 2));
        expect(!matrix.portrait().contains(0, 3));
        expect(!matrix.portrait().contains(1, 0));
        expect( matrix.portrait().contains(1, 1));
        expect(!matrix.portrait().contains(1, 2));
        expect( matrix.portrait().contains(1, 3));
        expect( matrix.portrait().contains(2, 0));
        expect(!matrix.portrait().contains(2, 1));
        expect(!matrix.portrait().contains(2, 2));
        expect(!matrix.portrait().contains(2, 3));
        expect(!matrix.portrait().contains(3, 0)); // Out of range row index
        expect(!matrix.portrait().contains(0, 4)); // Out of range column index

        // Check that shift() method returns correct indices for non-zero elements
        expect(eq(matrix.portrait().shift(0, 0), 0u));
        expect(eq(matrix.portrait().shift(0, 2), 1u));
        expect(eq(matrix.portrait().shift(1, 1), 2u));
        expect(eq(matrix.portrait().shift(1, 3), 3u));
        expect(eq(matrix.portrait().shift(2, 0), 4u));
        expect(throws<std::out_of_range>([&matrix] { matrix.portrait().shift(0, 1); })); // Non-zero element does not exist
        expect(throws<std::out_of_range>([&matrix] { matrix.portrait().shift(3, 0); })); // Out of range row index
        expect(throws<std::out_of_range>([&matrix] { matrix.portrait().shift(0, 4); })); // Out of range column index

        // Fill values and check that they are correctly accessed through operator()
        expect(nothrow([&matrix] { matrix(0, 0) = 1.0; }));
        expect(nothrow([&matrix] { matrix(0, 2) = 3.0; }));
        expect(nothrow([&matrix] { matrix(1, 1) = 2.0; }));
        expect(nothrow([&matrix] { matrix(1, 3) = 4.0; }));
        expect(nothrow([&matrix] { matrix(2, 0) = 5.0; }));

        expect(eq(matrix(0, 0), 1.0));
        expect(eq(matrix(0, 2), 3.0));
        expect(eq(matrix(1, 1), 2.0));
        expect(eq(matrix(1, 3), 4.0));
        expect(eq(matrix(2, 0), 5.0));
    };

    "validation"_test = [] {
        sparse_matrix<T, I, I> matrix{3, 3};

        // validate empty portrait
        expect(nothrow([&matrix] { validate_shifts(matrix.portrait().shifts()); }));

        // Invalid shifts size is less than 2
        matrix.portrait().shifts() = {0};
        expect(throws<std::logic_error>([&matrix] { validate_shifts(matrix.portrait().shifts()); }));

        // Invalid shifts start with non-zero
        matrix.portrait().shifts() = {1, 2, 3};
        expect(throws<std::logic_error>([&matrix] { validate_shifts(matrix.portrait().shifts()); }));

        // Invalid shifts vector contains negative values
        matrix.portrait().shifts() = {0, 2, -1, 3};
        expect(throws<std::logic_error>([&matrix] { validate_shifts(matrix.portrait().shifts()); }));

        // Invalid shifts vector is not non-decreasing
        matrix.portrait().shifts() = {0, 2, 1, 3};
        expect(throws<std::logic_error>([&matrix] { validate_shifts(matrix.portrait().shifts()); }));
        
        // Invalid indices size does not match the last element of shifts
        matrix.portrait().shifts() = {0, 2, 4, 6};
        matrix.portrait().indices() = {0, 2, 1, 2, 0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid shift size in a row is greater than the number of columns
        matrix.portrait().shifts() = {0, 4, 4, 5};
        matrix.portrait().indices() = {0, 1, 2, 3, 0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid column indices are negative
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, -1, 1, 2, 0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid column indices are out of range
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 1, 2, 3};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid column indices are not unique within a row
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 1, 1, 0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid column indices are not sorted within a row
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 2, 1, 0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix_portrait(matrix.portrait()); }));

        // Invalid values size does not match the number of non-zero elements in the portrait
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 1, 2, 0};
        matrix.values() = {1.0, 3.0, 2.0, 4.0};
        expect(throws<std::logic_error>([&matrix] { validate_sparse_matrix(matrix); }));

        // Valid matrix
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 1, 2, 0};
        matrix.values() = {1.0, 3.0, 2.0, 4.0, 5.0};
        expect(nothrow([&matrix] { validate_sparse_matrix(matrix); }));
    };

    "operations"_test = [] {
        static constexpr T Epsilon = std::numeric_limits<T>::epsilon();

        // Test matrix-vector multiplication for the following matrix and vector:
        // [1 0 3]   [1]   [10]
        // [0 2 4] * [2] = [16]
        // [5 0 0]   [3]   [ 5]
        sparse_matrix<T> matrix{3, 3};
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 1, 2, 0};
        matrix.values() = {1.0, 3.0, 2.0, 4.0, 5.0};
        const std::vector<T> vector = {1.0, 2.0, 3.0};
        const auto result = matrix * vector;
        expect(eq(result.size(), 3));
        expect(approx(result[0], 10.0, Epsilon));
        expect(approx(result[1], 16.0, Epsilon));
        expect(approx(result[2], 5.0, Epsilon));

        // Test matrix-vector multiplication for the following matrix and vector:
        // [1 0 3]   [1]   [10]
        // [0 2 4] * [2] = [16]
        // [3 4 0]   [3]   [11]
        const auto upper_triangular_result = matrix.self_adjoint<matrix_part::Upper>() * vector;
        expect(eq(upper_triangular_result.size(), 3));
        expect(approx(upper_triangular_result[0], 10.0, Epsilon));
        expect(approx(upper_triangular_result[1], 16.0, Epsilon));
        expect(approx(upper_triangular_result[2], 11.0, Epsilon));

        // Test matrix-vector multiplication for the following matrix and vector:
        // [1 0 5]   [1]   [16]
        // [0 2 0] * [2] = [ 4]
        // [5 0 0]   [3]   [ 5]
        const auto lower_triangular_result = matrix.self_adjoint<matrix_part::Lower>() * vector;
        expect(eq(lower_triangular_result.size(), 3));
        expect(approx(lower_triangular_result[0], 16.0, Epsilon));
        expect(approx(lower_triangular_result[1], 4.0, Epsilon));
        expect(approx(lower_triangular_result[2], 5.0, Epsilon));

        // Test block matrix-vector multiplication for the following block matrix and block vector:
        // [1 2|1 4|0 0]   [1]   [24]
        // [2 2|0 3|0 0]   [2]   [18]
        // [-----------]
        // [2 0|0 0|3 0] * [3] = [17]
        // [0 5|0 0|2 3]   [4]   [38]
        // [-----------]
        // [0 0|4 5|0 0]   [5]   [32]
        // [0 0|6 7|0 0]   [6]   [46]
        sparse_matrix<square_matrix<T, 2>> block_matrix{3, 3};
        block_matrix.portrait().shifts() = {0, 2, 4, 5};
        block_matrix.portrait().indices() = {0, 1, 0, 2, 1};
        block_matrix.values() = {{1.0, 2.0, 2.0, 2.0}, 
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
};

}