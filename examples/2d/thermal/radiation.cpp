#include "influence_functions_2d.hpp"
#include "nonstationary_heat_equation_solver_2d.hpp"
#include "heat_equation_solution_2d.hpp"

#include <iostream>

namespace {

using T = double;
using I = int64_t;
using namespace std::numbers;
using namespace metamath::symbolic;

// std::array<T, 2> Idx(const std::array<T, 2>& point, const T r) {
//     static constexpr T L = 1;
//     static const auto [x, y, lx, ly, a] = make_variables<5>();
//     static constexpr auto Ix = 
//         -exp(-(lx*lx + 2*y*y) / (2*a*a)) / (4*pi_v<T>) *
//         (erf<T>((lx - x) / (sqrt2_v<T> * a)) + erf<T>(x / (sqrt2_v<T> * a))) * 
//         (a * std::sqrt(2 * pi_v<T>) * (exp((y*(2*lx+y)) / (2*a*a)) - exp((lx*lx + y*y) / (2*a*a))) +
//         pi_v<T> * y * exp((lx*lx + 2*y*y) / (2*a*a)) * (-erf<T>(y/(sqrt2_v<T> * a)) + erf<T>((-lx+y)/(sqrt2_v<T> * a))));
//     static constexpr auto Ixdx = Ix.derivative<x>();
//     return {
//         Ixdx(point[0], point[1], L, L, r),
//         Ixdx(point[1], point[0], L, L, r),
//     };
// }

std::array<T, 2> Idx(const std::array<T, 2>& point) {
    static constexpr T L = 1;
    static const auto [x, y] = make_variables<2>();
    static constexpr auto Ix = (-6*y*y + 8*y - 6*x*x + 6*x + 19) / (192 * pi_v<T>);
    static constexpr auto Iy = (-6*x*x + 8*x - 6*y*y + 6*y + 19) / (192 * pi_v<T>);
    static constexpr auto Ixdx = Ix.derivative<x>();
    static constexpr auto Iydy = Iy.derivative<y>();
    return { Ixdx(point), Iydy(point) };
}

void save(const auto& mesh, const auto& parameters, const auto& temperature, const std::string folder, const size_t step) {
    auto solution = nonlocal::thermal::heat_equation_solution_2d{mesh, parameters, temperature};
    //solution.calc_flux();
    using namespace std::literals;
    //const auto& [TX, TY] = solution.flux();
    //solution.save_as_vtk(folder + "/"s + std::to_string(step) + "heat.vtk"s);
    nonlocal::mesh::utils::save_as_csv(folder + "/"s + std::to_string(step) + "T.csv"s, mesh->container(), temperature);
    //nonlocal::mesh::utils::save_as_csv(folder + "/"s + std::to_string(step) + "TX.csv"s, mesh->container(), TX);
    //nonlocal::mesh::utils::save_as_csv(folder + "/"s + std::to_string(step) + "TY.csv"s, mesh->container(), TY);
    std::cout << "step = " << step << std::endl;
}

}

int main(const int argc, const char *const *const argv) {
    static constexpr T p1 = 0.5;
    static constexpr T r = 2.;
    static const nonlocal::influence::normal_distribution_2d<T> influence{r};

    static constexpr T t0 = 100;
    static constexpr T L = 1;
    static constexpr T Tc = 1;
    static constexpr T kappa = 1;
    static constexpr T lambda = 1;
    static constexpr T alpha = 0;
    static constexpr T eps_r = 0.7;
    static constexpr T Ar = 0.7;

    auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
    mesh->find_neighbours({
        {"DEFAULT", 3 * r}
    });

    nonlocal::thermal::parameters_2d<T> parameters;
    parameters["DEFAULT"] = {
        .model = {
            .influence = influence,
            .local_weight = p1
        },
        .physical = {
            .conductivity = {lambda}
        }
    };

    static constexpr auto right_part = 
        [factor = Tc * lambda / (L*L), p2 = nonlocal::nonlocal_weight(p1)](const std::array<T, 2>& x) {
            const std::array<T, 2> dIdx = Idx(x);
            return factor * (kappa - p2 * (dIdx[0] + dIdx[1]));
        };

    static constexpr auto init_dist = [](const std::array<T, 2>& x) constexpr noexcept {
        return x[0] * x[1] + t0;
    };

    nonlocal::thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::combined_flux_2d<T>>(
        [p2 = nonlocal::nonlocal_weight(p1)](const std::array<T, 2>& x) {
            return -p1 * x[1] - p2 * Idx(x)[0] - alpha * (Tc - t0);
        },
        alpha, Tc,
        eps_r
    );
    boundaries_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        [p2 = nonlocal::nonlocal_weight(p1)](const std::array<T, 2>& x) {
            return p1 * x[1] + p2 * Idx(x)[0];
        });
    boundaries_conditions["Up"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        [p2 = nonlocal::nonlocal_weight(p1)](const std::array<T, 2>& x) {
            return p1 * x[0] + p2 * Idx(x)[1];
        });
    boundaries_conditions["Down"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        [p2 = nonlocal::nonlocal_weight(p1)](const std::array<T, 2>& x) {
            return -p1 * x[0] - p2 * Idx(x)[1];
        });

    static constexpr T tau = 0.001;
    nonlocal::thermal::nonstationary_heat_equation_solver_2d<T, I, I> solver{mesh, tau};
    solver.compute(parameters, boundaries_conditions, init_dist);
    save(mesh, parameters, solver.temperature(), argv[2], 0u);
    for(const size_t step : std::ranges::iota_view{1u, 1001u}) {
        boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::combined_flux_2d<T>>(
            [p2 = nonlocal::nonlocal_weight(p1), t = t0 + tau * step](const std::array<T, 2>& x) {
                using nonlocal::thermal::STEFAN_BOLTZMANN_CONSTANT;
                using metamath::functions::power;
                return -p1 * x[1] - p2 * Idx(x)[0] - alpha * (Tc - t) + eps_r * STEFAN_BOLTZMANN_CONSTANT<T> * power<4>(t);
            },
            alpha, Tc,
            eps_r
        );
        solver.calc_step(boundaries_conditions, right_part);
        if (step % 500 == 0)
            save(mesh, parameters, solver.temperature(), argv[2], step);
    }

    return 0;
}