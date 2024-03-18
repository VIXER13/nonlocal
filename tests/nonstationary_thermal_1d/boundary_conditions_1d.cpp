#include "nonstationary_tests_1d.hpp"

namespace bc_1d_tests {

using namespace nonstat_1d_tests;
using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

const suite<"boundary_conditions_temperature_flux"> _ = [] {
    "Test №1"_test = [] {
        // Two first kind BCs
        // BC T|x=0 = alpha, T|x=L = beta
        // Exact solution T(x) = alpha * (x - L) / L + beta * x / L
        std::function<T(T)> ref_sol = [](T x) {
            return 10.0 * (2.0 - x) / 2.0 - 3.0 * x / 2.0;
        };
        std::filesystem::path path = "/mnt/c/repos/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test1.json";
        solve_nonstationary_thermal_1d_problem<T, I>(nonlocal::config::parse_json(path), ref_sol);
    };

    "Test №2"_test = [] {
        // Zero flux with temperature
        // T|x=0 = alpha, q*n|x=L = 0 OR q*n|x=0 = 0, T|x=L = alpha
        // Exact solution T(x) = alpha
        std::function<T(T)> ref_sol = [](T x) {
            return 10.0;
        };
        std::filesystem::path path1 = "/mnt/c/repos/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test2_1.json";
        std::filesystem::path path2 = "/mnt/c/repos/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test2_2.json";
        solve_nonstationary_thermal_1d_problem<T, I>(nonlocal::config::parse_json(path1), ref_sol);
        solve_nonstationary_thermal_1d_problem<T, I>(nonlocal::config::parse_json(path2), ref_sol);
    };
};
    
}
