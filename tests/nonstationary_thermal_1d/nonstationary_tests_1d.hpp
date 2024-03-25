#ifndef UNIT_TESTS_BOUNDARY_CONDITIONS_HPP
#define UNIT_TESTS_BOUNDARY_CONDITIONS_HPP

#include <boost/ut.hpp>
#include <tuple>
#include <filesystem>

#include "problems_utils.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

namespace nonstat_1d_tests {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

template <std::floating_point T>
struct time_data final {
    T time_step = T{0};       
    T initial_time = T{0};
    uint64_t steps_count = 0; 
    uint64_t save_frequency = 1ull;

    explicit constexpr time_data() noexcept = default;
    explicit time_data(T _time_step, T _initial_time, uint64_t _steps_count, uint64_t _save_frequency)
        : time_step(_time_step), initial_time(_initial_time), steps_count(_steps_count), save_frequency(_save_frequency) {};
};

template<std::floating_point T, std::signed_integral I>
void check_solution(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const heat_equation_solution_1d<T>& solution, T time_layer, I step, 
                    std::function<T(T, T)> ref, T eps = epsilon) {
    auto& sol = solution.temperature();
    for (std::size_t k = 0; k < sol.size(); ++k) {
        expect(lt(std::fabs(sol[k] - ref(time_layer, mesh->node_coord(k))), eps * step));
        std::cout << "step: " <<  step << std::endl;
    }
}

template<std::floating_point T, std::signed_integral I>
void save_and_calc_flux(const nonstat_1d_tests::time_data<T>& time, I step, heat_equation_solution_1d<T>& solution) {
    if (step % time.save_frequency == 0) {
        solution.calc_flux();
        const std::filesystem::path path = "/results";
        mesh::utils::save_as_csv(path, solution.mesh(), {{"temperature", solution.temperature()}, {"flux", solution.flux()}});
    }
}
    
}

#endif
