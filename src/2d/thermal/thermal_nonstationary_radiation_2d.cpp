#include "nonstationary_heat_equation_solver_2d.hpp"
#include "heat_equation_solution_2d.hpp"
#include "influence_functions_2d.hpp"
#include "mesh_container_2d_utils.hpp"

#include "metamath.hpp"

#include <iostream>
#include <cmath>

namespace {

using T = double;
using I = int64_t;

void save_solution(const nonlocal::thermal::heat_equation_solution_2d<T, I>& solution, 
                   const std::filesystem::path& folder, const size_t step) {
    solution.save_as_vtk(folder.string() + '/' + std::to_string(step) + "T.vtk");
    nonlocal::mesh::utils::save_as_csv(folder.string() + '/' + std::to_string(step) + "T.csv", solution.mesh().container(), solution.temperature());
    if (solution.is_flux_calculated()) {
        const auto& [TX, TY] = solution.flux();
        nonlocal::mesh::utils::save_as_csv(folder.string() + '/' + std::to_string(step) + "TX.csv", solution.mesh().container(), TX);
        nonlocal::mesh::utils::save_as_csv(folder.string() + '/' + std::to_string(step) + "TY.csv", solution.mesh().container(), TY);
    }
}

}

int main(int argc, char** argv) {
#ifdef MPI_USE
    MPI_Init(&argc, &argv);
#endif

    if(argc < 6) {
        std::cerr << "Input format [program name] <path to mesh> <r1> <r2> <p1> <save_path>" << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);

        const std::array<T, 2> r = {std::stod(argv[2]), std::stod(argv[3])};

        static const nonlocal::influence::normal_distribution_2d<T> influence_function{r};

        const T p1 = std::stod(argv[4]);

        static constexpr T tau = 0.001;
        nonlocal::thermal::parameter_2d<T> parameters = {
            .conductivity = {T{1}, T{0},
                             T{0}, T{1}},
            .capacity = T{1},
            .density = T{1},
            .material = nonlocal::material_t::ISOTROPIC
        };

        if (nonlocal::theory_type(p1) == nonlocal::theory_t::NONLOCAL) {
            mesh->find_neighbours(std::max(r[0], r[1]));
        };

        const T pi = M_PI;
        const T lambda = T{1};
        constexpr T ambient_temperature = 100;
        const T sigma = nonlocal::thermal::STEFAN_BOLTZMANN_CONSTANT<T>;
        const T p2 = 1 - p1;
        const std::vector<T> coefs = {-0.5 / (sqrt(2 * pi)), 
                                       1 / (sqrt(2) * r[0]), 
                                       0.5 / (r[0] * r[0]),
                                       sqrt(pi * 0.5),
                                       0.25 / sqrt(r[0] * pi),
                                       2 / sqrt(pi)};

        auto twoexps = [](T x, T coef)
        {
            return std::exp(x * x * coef) - std::exp((x - 1) * (x - 1) * coef);
        };

        auto twoerfs = [](T x, T coef)
        {
            return std::erf(x * coef) - std::erf((x - 1) * coef);
        };

        auto Ifunc = [&r, &coefs, twoexps, twoerfs](const std::array<T, 2>& x) constexpr noexcept { 
            return coefs[0] * twoerfs(x[0], coefs[1]) * 
            (r[0] * twoexps(x[1], coefs[2]) + x[1] * coefs[3] * twoerfs(x[1], coefs[1]) ); 
        };

        auto divI = [&r, &coefs, twoexps, twoerfs](const std::array<T, 2>& x)
        {
            return coefs[4] * (
                   coefs[5] * twoexps(x[0], coefs[2]) * 
                   (r[0] * twoexps(x[1], coefs[2]) + x[1] * coefs[3] * twoerfs(x[1], coefs[1])) +
                   coefs[5] * twoexps(x[1], coefs[2]) * 
                   (r[0] * twoexps(x[0], coefs[2]) + x[0] * coefs[3] * twoerfs(x[0], coefs[1]))
                   );
        };

        const std::vector<T> absorption = {0.7, 0.0, 0.0, 0.0}; //Left, Right, Up, Down
        const std::vector<T> emissivity = {0.7, 0.0, 0.0, 0.0};
        const std::vector<T> heat_transfer = {0.7, 0.0, 0.0, 0.0};
        const T lambda_temp_Ar = lambda * ambient_temperature / absorption[0];
        const T Bi = heat_transfer[0] / lambda;
        const T N_er =  emissivity[0] * sigma * metamath::functions::power<3>(ambient_temperature) / lambda;

        // auto left_q = [lambda_temp_Ar, N_er, Bi, p1, p2, Ifunc, t](const std::array<T, 2>& x)
        // {
        //     return lambda_temp_Ar * (N_er * metamath::functions::power<4>(t) - Bi * (1 - t) 
        //            - p1 * x[1] - p2 * Ifunc({x[0], x[1]}));
        // };

        auto left_q = [](const std::array<T, 2>& x)
        {
            return T{1};
        };

        nonlocal::thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
        boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::combined_flux_2d<T>>(
            left_q, heat_transfer[0], ambient_temperature, emissivity[0]
        );
        boundaries_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
            [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
        );
        boundaries_conditions["Up"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
            [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
        );
        boundaries_conditions["Down"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
            [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
        );

        static constexpr auto init_dist = [](const std::array<T, 2>& x) constexpr noexcept {
            return x[0] * x[1];
        };

        const T c_rho_temp = T{1} * T{1} * ambient_temperature; //c * rho * Tc
        static auto right_part = [c_rho_temp, lambda, p2, divI, ambient_temperature](const std::array<T, 2>& x) constexpr noexcept {
            return c_rho_temp - ambient_temperature * lambda * p2 * divI({x[0], x[1]}) ;
        };
        
        const std::filesystem::path FOLDER =  argv[5];

        nonlocal::thermal::nonstationary_heat_equation_solver_2d<T, I, I> solver{mesh, tau};
        solver.compute(parameters, boundaries_conditions, init_dist, p1, influence_function);
        for(const size_t step : std::ranges::iota_view{0u, 101u}) {
            solver.calc_step(boundaries_conditions, right_part);
            auto solution = nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, p1, influence_function, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, FOLDER, step);
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

#ifdef MPI_USE
    MPI_Finalize();
#endif
    return 0;
}