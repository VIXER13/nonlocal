#include <solvers/solver_2d/mechanical/elastic_parameters.hpp>

#include <boost/ut.hpp>

#include <numbers>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal;
using namespace solver_2d::mechanical;
using namespace std::numbers;

const suite<"elastic_parameters_2d"> _ = [] {
    "isotropic"_test = [] {
        const orthotropic_elastic_parameters<T> elastic = {
            .young_modulus = {350., 350.},
            .poissons_ratio = {0.25, 0.25},
            .shear_modulus = 140.
        };
        const auto hooke = elastic.hooke({});
        for(const T angle : {T{0}, pi / 6, -pi / 6, pi / 4, -pi / 4, pi / 2, pi}) {
            using namespace anisotropic_indices;
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            static constexpr T Epsilon = std::numeric_limits<T>::epsilon();
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
        static constexpr T Epsilon = 2.5e-15;
        for(const T angle : {T{0}, pi}) {
            using namespace anisotropic_indices;
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            expect(approx(result[_11], hooke[_11], Epsilon));
            expect(approx(result[_12], hooke[_12], Epsilon));
            expect(approx(result[_16],       T{0}, Epsilon));
            expect(approx(result[_22], hooke[_22], Epsilon));
            expect(approx(result[_26],       T{0}, Epsilon));
            expect(approx(result[_66], hooke[_66], Epsilon));
        }

        for(const T angle : {pi / 2, -pi / 2}) {
            using namespace anisotropic_indices;
            const auto result = anisotropic_elastic_parameters<T>::rotate(hooke, angle);
            expect(approx(result[_11], hooke[_22], Epsilon));
            expect(approx(result[_12], hooke[_12], Epsilon));
            expect(approx(result[_16],       T{0}, Epsilon));
            expect(approx(result[_22], hooke[_11], Epsilon));
            expect(approx(result[_26],       T{0}, Epsilon));
            expect(approx(result[_66], hooke[_66], Epsilon));
        }
    };
};

}