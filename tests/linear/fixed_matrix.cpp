#include <metamath/linear/fixed_matrix.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;

using T = double;

constexpr T Epsilon = std::numeric_limits<T>::epsilon();

suite<"fixed_matrix"> _ = [] {
    "determinant"_test = [] {
        static constexpr square_matrix<T, 1> Matrix_1x1 = {5.0};
        expect(eq(determinant(Matrix_1x1), 5.0));

        static constexpr square_matrix<T, 2> Matrix_2x2 = {1.0, 2.0, 
                                                           3.0, 4.0};
        expect(approx(determinant(Matrix_2x2), -2.0, Epsilon));

        static constexpr square_matrix<T, 3> Matrix_3x3 = {6.0,  1.0, 1.0, 
                                                           4.0, -2.0, 5.0, 
                                                           2.0,  8.0, 7.0};
        expect(approx(determinant(Matrix_3x3), -306.0, Epsilon));
    };

    "is_positive"_test = [] {
        static constexpr square_matrix<T, 1> Matrix_1x1_Positive = {5.0};
        expect(is_positive(Matrix_1x1_Positive));

        static constexpr square_matrix<T, 1> Matrix_1x1_Negative = {-5.0};
        expect(!is_positive(Matrix_1x1_Negative));

        static constexpr square_matrix<T, 1> Matrix_1x1_Zero = {0.0};
        expect(!is_positive(Matrix_1x1_Zero));

        static constexpr square_matrix<T, 2> Matrix_2x2_Positive = {4.0, 1.0, 
                                                                    1.0, 3.0};
        expect(is_positive(Matrix_2x2_Positive));

        static constexpr square_matrix<T, 2> Matrix_2x2_Negative = {1.0, 2.0, 
                                                                    3.0, 4.0};
        expect(!is_positive(Matrix_2x2_Negative));

        static constexpr square_matrix<T, 3> Matrix_3x3_Positive = {6.0, 1.0, 1.0, 
                                                                    1.0, 4.0, 5.0, 
                                                                    1.0, 5.0, 7.0};
        expect(is_positive(Matrix_3x3_Positive));

        static constexpr square_matrix<T, 3> Matrix_3x3_Negative = {6.0,  1.0, 1.0, 
                                                                    4.0, -2.0, 5.0, 
                                                                    2.0,  8.0, 7.0};
        expect(!is_positive(Matrix_3x3_Negative));
    };

    "inverse"_test = [] {
        static constexpr square_matrix<T, 1> Matrix_1x1 = {5.0};
        static constexpr square_matrix<T, 1> Expected_1x1_Inverse = {0.2};
        const auto inverse_1x1 = inverse(Matrix_1x1);
        expect(approx(inverse_1x1[0][0], Expected_1x1_Inverse[0][0], Epsilon));

        static constexpr square_matrix<T, 2> Matrix_2x2 = {4.0, 1.0, 
                                                           1.0, 3.0};
        static constexpr square_matrix<T, 2> Expected_2x2_Inverse = { 3.0 / 11.0, -1.0 / 11.0, 
                                                                     -1.0 / 11.0,  4.0 / 11.0};
        static constexpr auto Inverse_2x2 = inverse(Matrix_2x2);
        for (const size_t row : std::ranges::iota_view{0zu, 2zu})
            for (const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(Inverse_2x2[row][col], Expected_2x2_Inverse[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";

        static constexpr square_matrix<T, 3> Matrix_3x3 = {6.0,  1.0, 1.0, 
                                                           4.0, -2.0, 5.0, 
                                                           2.0,  8.0, 7.0};
        static constexpr square_matrix<T, 3> Expected_3x3_Inverse = { 54.0 / 306.0,  -1.0 / 306.0, -7.0 / 306.0,
                                                                      18.0 / 306.0, -40.0 / 306.0, 26.0 / 306.0,
                                                                     -36.0 / 306.0,  46.0 / 306.0, 16.0 / 306.0};
        static constexpr auto Inverse_3x3 = inverse(Matrix_3x3);
        for (const size_t row : std::ranges::iota_view{0zu, 3zu})
            for (const size_t col : std::ranges::iota_view{0zu, 3zu})
                expect(approx(Inverse_3x3[row][col], Expected_3x3_Inverse[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "multiplication"_test = [] {
        static constexpr fixed_matrix<T, 2, 3> Matrix_2x3 = {1.0, 2.0, 3.0,
                                                             4.0, 5.0, 6.0};
        static constexpr std::array<T, 3> Vector_3 = {7.0, 8.0, 9.0};
        static constexpr std::array<T, 2> Expected_Product_2 = {50.0, 122.0};
        static constexpr auto Product_2 = Matrix_2x3 * Vector_3;
        for (const size_t i : std::ranges::iota_view{0zu, 2zu})
            expect(approx(Product_2[i], Expected_Product_2[i], Epsilon)) << " at element " << i;

        static constexpr fixed_matrix<T, 3, 2> Matrix_3x2 = {1.0, 4.0,
                                                             2.0, 5.0,
                                                             3.0, 6.0};
        static constexpr std::array<T, 2> Vector_2 = {7.0, 8.0};
        static constexpr std::array<T, 3> Expected_Product_3 = {39.0, 54.0, 69.0};
        static constexpr auto Product_3 = Matrix_3x2 * Vector_2;
        for (const size_t i : std::ranges::iota_view{0zu, 3zu})
            expect(approx(Product_3[i], Expected_Product_3[i], Epsilon)) << " at element " << i;

        static constexpr fixed_matrix<T, 2, 3> Matrix_A = {1.0, 2.0, 3.0,
                                                           4.0, 5.0, 6.0};
        static constexpr fixed_matrix<T, 3, 2> Matrix_B = { 7.0,  8.0,
                                                            9.0, 10.0,
                                                           11.0, 12.0};
        static constexpr fixed_matrix<T, 2, 2> Expected_Product = { 58.0,  64.0,
                                                                      139.0, 154.0};
        static constexpr auto Product_AB = Matrix_A * Matrix_B;
        for (const size_t row : std::ranges::iota_view{0zu, 2zu})
            for (const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(Product_AB[row][col], Expected_Product[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "multiplication_assignment"_test = [] {
        fixed_matrix<T, 2, 2> Matrix = {1.0, 2.0,
                                        3.0, 4.0};
        static constexpr fixed_matrix<T, 2, 2> Multiplier = {5.0, 6.0,
                                                             7.0, 8.0};
        static constexpr fixed_matrix<T, 2, 2> Expected_Product = {19.0, 22.0,
                                                                   43.0, 50.0};
        Matrix *= Multiplier;
        for (const size_t row : std::ranges::iota_view{0zu, 2zu})
            for (const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(Matrix[row][col], Expected_Product[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "matrix_division"_test = [] {
        static constexpr fixed_matrix<T, 2, 2> Matrix = {2.0, 4.0,
                                                         6.0, 8.0};
        static constexpr fixed_matrix<T, 2, 2> Divisor = {1.0, 2.0,
                                                          3.0, 4.0};
        static constexpr fixed_matrix<T, 2, 2> Expected_Quotient = {2.0, 0.0,
                                                                    0.0, 2.0};
        const auto Quotient = Matrix / Divisor;
        for (const size_t row : std::ranges::iota_view{0zu, 2zu})
            for (const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(Quotient[row][col], Expected_Quotient[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };

    "division_assignment"_test = [] {
        fixed_matrix<T, 2, 2> Matrix = {2.0, 4.0,
                                        6.0, 8.0};
        static constexpr fixed_matrix<T, 2, 2> Divisor = {1.0, 2.0,
                                                          3.0, 4.0};
        static constexpr fixed_matrix<T, 2, 2> Expected_Quotient = {2.0, 0.0,
                                                                    0.0, 2.0};
        Matrix /= Divisor;
        for (const size_t row : std::ranges::iota_view{0zu, 2zu})
            for (const size_t col : std::ranges::iota_view{0zu, 2zu})
                expect(approx(Matrix[row][col], Expected_Quotient[row][col], Epsilon)) << " at element (" << row << ", " << col << ")";
    };
};

}