#include <metamath/linear/sparce_matrix.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;

using T = double;

sparce_matrix<T> init_matrix() {
    // Matrix representation:
    // [1 0 3]
    // [0 2 4]
    // [5 0 0]
    sparce_matrix<T> matrix{3, 3};
    matrix.portrait().shifts() = {0, 2, 4, 5};
    matrix.portrait().indices() = {0, 2, 1, 2, 0};
    matrix.values() = {1.0, 3.0, 2.0, 4.0, 5.0};
    return matrix;
}

suite<"sparce_matrix"> _ = [] {
    "common"_test = [] {
        const auto matrix = init_matrix();

        expect(matrix.rows() == 3);
        expect(matrix.cols() == 3);
        expect(matrix.non_zeros() == 5);

        expect( matrix.portrait().contains(0, 0));
        expect(!matrix.portrait().contains(0, 1));
        expect( matrix.portrait().contains(0, 2));
        expect(!matrix.portrait().contains(1, 0));
        expect( matrix.portrait().contains(1, 1));
        expect( matrix.portrait().contains(1, 2));
        expect( matrix.portrait().contains(2, 0));
        expect(!matrix.portrait().contains(2, 1));
        expect(!matrix.portrait().contains(2, 2));

        expect(eq(matrix.portrait().shift(0, 0), 0u));
        expect(eq(matrix.portrait().shift(0, 2), 1u));
        expect(eq(matrix.portrait().shift(1, 1), 2u));
        expect(eq(matrix.portrait().shift(1, 2), 3u));
        expect(eq(matrix.portrait().shift(2, 0), 4u));

        expect(eq(matrix(0, 0), 1.0));
        expect(eq(matrix(0, 2), 3.0));
        expect(eq(matrix(1, 1), 2.0));
        expect(eq(matrix(1, 2), 4.0));
        expect(eq(matrix(2, 0), 5.0));
    };

    "validation"_test = [] {
        sparce_matrix<T> matrix{3, 3};

        // Valid empty portrait
        matrix.portrait().shifts() = {};
        expect(nothrow([&matrix] { matrix.portrait().validate(); }));

        // Invalid shifts size is less than 2
        matrix.portrait().shifts() = {0};
        expect(throws<std::logic_error>([&matrix] { matrix.portrait().validate(); }));
        
        // Invalid indices size does not match the last element of shifts
        matrix.portrait().shifts() = {0, 2, 4, 6};
        matrix.portrait().indices() = {0, 2, 1, 2, 0};
        expect(throws<std::logic_error>([&matrix] { matrix.portrait().validate(); }));

        // Invalid shifts vector is not non-decreasing
        matrix.portrait().shifts() = {0, 4, 2, 5};
        expect(throws<std::logic_error>([&matrix] { matrix.portrait().validate(); }));

        // Invalid column indices are out of range
        matrix.portrait().shifts() = {0, 2, 4, 6};
        matrix.portrait().indices() = {0, 2, 1, 2, 3};
        expect(throws<std::logic_error>([&matrix] { matrix.portrait().validate(); }));

        // Invalid column indices are not sorted within a row
        matrix.portrait().shifts() = {0, 2, 4, 5};
        matrix.portrait().indices() = {0, 2, 2, 1, 0};
        expect(throws<std::logic_error>([&matrix] { matrix.portrait().validate(); }));

        // Valid matrix
        matrix = init_matrix();
        expect(nothrow([&matrix] { matrix.portrait().validate(); }));
        expect(nothrow([&matrix] { matrix.validate(); }));
    };
};

}