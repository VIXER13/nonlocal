#include "nonstationary_heat_equation_solver_2d.hpp"
#include "heat_equation_solution_2d.hpp"
#include "influence_functions_2d.hpp"
#include "mesh_container_2d_utils.hpp"

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

template<class T>
T Ifunc(const T r, const std::array<T, 2>& x) noexcept {
    static constexpr T pi = std::numbers::pi_v<double>;
    static constexpr T sqrt2_v = std::numbers::sqrt2_v<double>;
    return T{1} / (4 * pi) * std::exp(-(1 + 2 * x[0]*x[0]) / (2 * r * r)) *
           (std::erf((x[1] - 1) / (sqrt2_v * r)) - std::erf(x[1] / sqrt2_v)) *
           (
                std::sqrt(2 * pi) * r * (std::exp(x[0] * (2 + x[0]) / (2 * r * r)) - std::exp((1 + x[0] * x[0])) / (2 * r * r)) +
                pi * x[0] * std::exp((1 + 2 * x[0]*x[0]) / (2 * r * r)) *
                    (std::erf((x[0] - 1) / (sqrt2_v * r)) - std::erf(x[0] / (sqrt2_v * r)))
           );
}

template<class T>
T flux_x(const std::array<T, 2>& p, const T r, const std::array<T, 2>& x) noexcept {
    return -p[0] * x[1];// + p[1] * Ifunc(r, {x[0], x[1]});
}

template<class T>
T flux_y(const std::array<T, 2>& p, const T r, const std::array<T, 2>& x) noexcept {
    return -p[0] * x[0];// + p[1] * Ifunc(r, {x[1], x[0]});
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

        static constexpr T tau = 0.01;
        nonlocal::thermal::parameter_2d<T> parameters = {
            .conductivity = {T{1}, T{0},
                             T{0}, T{1}},
            .capacity = T{1},
            .density = T{1},
            .material = nonlocal::material_t::ISOTROPIC
        };

        if (nonlocal::theory_type(p1) == nonlocal::theory_t::NONLOCAL) {
            mesh->find_neighbours(3 * std::max(r[0], r[1]));
        }

        nonlocal::thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
        boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(1, 5);
        boundaries_conditions["Right"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(1, 5);
        boundaries_conditions["Up"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(1, 5);
        boundaries_conditions["Down"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(1, 5);

        // boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [p1, rr = std::max(r[0], r[1])](const std::array<T, 2>& x) constexpr noexcept { return -flux_x({p1, 1-p1}, rr, x); }
        // );
        // boundaries_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [p1, rr = std::max(r[0], r[1])](const std::array<T, 2>& x) constexpr noexcept { return  flux_x({p1, 1-p1}, rr, x); }
        // );
        // boundaries_conditions["Up"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [p1, rr = std::max(r[0], r[1])](const std::array<T, 2>& x) constexpr noexcept { return  flux_y({p1, 1-p1}, rr, x); }
        // );
        // boundaries_conditions["Down"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [p1, rr = std::max(r[0], r[1])](const std::array<T, 2>& x) constexpr noexcept { return -flux_y({p1, 1-p1}, rr, x); }
        // );

        static constexpr auto init_dist = [](const std::array<T, 2>& x) constexpr noexcept {
            return 0;
        };

        static constexpr auto right_part = [](const std::array<T, 2>& x) constexpr noexcept {
            return 0;
        };
        
        nonlocal::thermal::nonstationary_heat_equation_solver_2d<T, I, I> solver{mesh, tau};
        solver.compute(parameters, boundaries_conditions, init_dist, p1, influence_function);
        for(const size_t step : std::ranges::iota_view{1u, 1001u}) {
            solver.calc_step(boundaries_conditions, right_part);
            auto solution = nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, p1, influence_function, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, "./sol", step);
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