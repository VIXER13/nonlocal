#include <config/read_mechanical_boundary_conditions.hpp>
#include <config/read_thermal_boundary_conditions.hpp>
#include <solvers/solver_1d/base/boundary_conditions_1d.hpp>

#include <boost/ut.hpp>

#include <array>
#include <ranges>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal::config;
using namespace nonlocal::solver_1d;

void check_robin_coefficients(const boundaries_data_1d<T>& data) {
    static constexpr auto Expected_Coefficients = std::array{ T{ 3 }, T{ 4 } };
    static constexpr auto Expected_Loads = std::array{ T{ 30 }, T{ 22 } };
    for (const size_t boundary : std::ranges::iota_view{ 0zu, 2zu }) {
        expect(data[boundary].robin_coefficient.has_value()) << fatal;
        expect(eq(*data[boundary].robin_coefficient, Expected_Coefficients[boundary]));
        expect(eq(data[boundary].load, Expected_Loads[boundary]));
    }
}

const suite<"config_boundary_conditions_1d"> _ = [] {
    "thermal_robin_coefficients"_test = [] {
        const nlohmann::json thermal_config = {
            { "left", { { "kind", "convection" }, { "temperature", 10. }, { "heat_transfer", 3. } } },
            { "right", { { "kind", "combined" }, { "flux", 2. }, { "temperature", 5. },
                           { "heat_transfer", 4. }, { "emissivity", 0.5 } } }
        };
        // alpha = h = 3, beta = h*T_ambient = 30; on the right, beta = 2 + 4*5 = 22.
        check_robin_coefficients(normalize_boundaries_1d(read_thermal_boundaries_conditions_1d<T>(thermal_config, {})));
    };

    "mechanical_robin_coefficients"_test = [] {
        const nlohmann::json mechanical_config = {
            { "left", { { "kind", "spring" }, { "displacement", 10. }, { "stiffness", 3. } } },
            { "right", { { "kind", "combined" }, { "pressure", 2. }, { "displacement", 5. }, { "stiffness", 4. } } }
        };
        // alpha = k = 3, beta = k*u_ref = 30; on the right, beta = p + k*u_ref = 2 + 4*5 = 22.
        check_robin_coefficients(normalize_boundaries_1d(read_mechanical_boundaries_conditions_1d<T>(mechanical_config, {})));
    };
};

} // namespace
