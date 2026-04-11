#include <solvers/solver_2d/mechanical/elastic_parameters.hpp>

#include <boost/ut.hpp>

#include <numbers>

namespace {

using T = double;
using namespace nonlocal;
using namespace solver_2d::mechanical;
using namespace anisotropic_indices;
using namespace boost::ut;
using namespace std::numbers;

const suite<"elastic_parameters_2d"> _ = [] {
    "isotropic"_test = [] {
        static constexpr T Epsilon = std::numeric_limits<T>::epsilon();
        const orthotropic_elastic_parameters<T> elastic = {
            .young_modulus = {350., 350.},
            .poissons_ratio = {0.25, 0.25},
            .shear_modulus = 140.
        };
        const auto hooke = elastic.hooke({});
        expect(approx(hooke[_11], hooke[_22], Epsilon));
        for(const size_t i : std::ranges::iota_view{0zu, 16zu}) {
            const T angle = -pi + i * (2 * pi) / 16;
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            expect(approx(result[_11], hooke[_11], Epsilon));
            expect(approx(result[_12], hooke[_12], Epsilon));
            expect(approx(result[_16],       T{0}, Epsilon));
            expect(approx(result[_22], hooke[_22], Epsilon));
            expect(approx(result[_26],       T{0}, Epsilon));
            expect(approx(result[_66], hooke[_66], Epsilon));
        }
    };

    "orthotropic"_test = [] {
        const orthotropic_elastic_parameters<T> elastic = {
            .young_modulus = {35., 21.},
            .poissons_ratio = {0.25, 0.15},
            .shear_modulus = 18.
        };
        const auto hooke = elastic.hooke({});

        for(const T angle : {T{0}, pi, pi / 2, -pi / 2}) {
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            static constexpr T Epsilon = 2.5e-15;
            if (angle == 0 || angle == pi) {
                expect(approx(result[_11], hooke[_11], Epsilon));
                expect(approx(result[_22], hooke[_22], Epsilon));
            } else {
                expect(approx(result[_11], hooke[_22], Epsilon));
                expect(approx(result[_22], hooke[_11], Epsilon));
            }
            expect(approx(result[_12], hooke[_12], Epsilon));
            expect(approx(result[_16],       T{0}, Epsilon));
            expect(approx(result[_26],       T{0}, Epsilon));
            expect(approx(result[_66], hooke[_66], Epsilon));
        }

        for(const T angle : {pi / 4, -pi / 4, -3 * pi / 4, 3 * pi / 4}) {
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            static constexpr T Epsilon = 7.2e-15;
            expect(approx(result[_11], result[_22], Epsilon));
            expect(approx(result[_16], result[_26], Epsilon));
            expect(approx(result[_11], T{0.25} * (hooke[_11] + hooke[_22] + 2 * hooke[_12]) + hooke[_66], Epsilon));
            expect(approx(result[_12], T{0.25} * (hooke[_11] + hooke[_22] + 2 * hooke[_12]) - hooke[_66], Epsilon));
            expect(approx(result[_66], T{0.25} * (hooke[_11] + hooke[_22] - 2 * hooke[_12]), Epsilon));
            if (angle == pi / 4 || angle == -3 * pi / 4)
                expect(approx(result[_16], T{0.25} * (hooke[_22] - hooke[_11]), Epsilon));
            else
                expect(approx(result[_16], T{0.25} * (hooke[_11] - hooke[_22]), Epsilon));
        }

        for(const auto [left_angle, right_angle] : {
                std::pair{pi / 7,     -pi / 7}, std::pair{-3 * pi / 5,  3 * pi / 5}, // OX Symmetry
                std::pair{pi / 7,  6 * pi / 7}, std::pair{-2 * pi / 9, -7 * pi / 9}, // OY Symmetry
                std::pair{pi / 7, -6 * pi / 7}, std::pair{-2 * pi / 5,  3 * pi / 5}  // Non-adjacent sectors
            }) {
            static constexpr T Epsilon = 1.5e-14;
            const auto left_result = anisotropic_elastic_parameters<T>::rotate(hooke, left_angle);
            const auto right_result = anisotropic_elastic_parameters<T>::rotate(hooke, right_angle);
            expect(approx(left_result[_11], right_result[_11], Epsilon));
            expect(approx(left_result[_12], right_result[_12], Epsilon));
            expect(approx(left_result[_66], right_result[_66], Epsilon));
            expect(approx(left_result[_22], right_result[_22], Epsilon));
            if (approx(std::abs(left_angle - right_angle), pi, Epsilon)) { // Non-adjacent sectors
                expect(approx(left_result[_16], right_result[_16], Epsilon));
                expect(approx(left_result[_26], right_result[_26], Epsilon));
            } else {
                expect(approx(left_result[_16], -right_result[_16], Epsilon));
                expect(approx(left_result[_26], -right_result[_26], Epsilon));
            }
        }
    };
};

}