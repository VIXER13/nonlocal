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
        static const nonlocal::influence::polynomial_2d<T, 2, 1> influence_function{r};
        const T p1 = std::stod(argv[4]);

        static constexpr T tau = 0.001;
        nonlocal::thermal::parameters_2d<T> parameters;
        parameters["Material1"] = {
            .model = {
                .influence = nonlocal::influence::polynomial_2d<T, 2, 1>{r},
                .local_weight = T{1}
            },
            .physical = {
                .conductivity = {T{1}}
            }
        };
        parameters["Material2"] = {
            .model = {
                .influence = nonlocal::influence::normal_distribution_2d<T>{std::stod(argv[2])},
                .local_weight = p1
            },
            .physical = {
                .conductivity = {T{10}}
            }
        };
        parameters["Material3"] = {
            .model = {
                .influence = nonlocal::influence::polynomial_2d<T, 2, 1>{r},
                .local_weight = T{1}
            },
            .physical = {
                .conductivity = {T{2}}
            }
        };

        mesh->find_neighbours({
            {"Material2", 3 * std::stod(argv[2])}
        });

        nonlocal::thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
        boundaries_conditions["Left"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(-1);
        boundaries_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(1);

        static constexpr auto init_dist = [](const std::array<T, 2>& x) constexpr noexcept {
            return T{0};
        };

        static constexpr auto right_part = [](const std::array<T, 2>& x) constexpr noexcept {
            return T{0};
        };
        
        nonlocal::thermal::nonstationary_heat_equation_solver_2d<T, I, I> solver{mesh, tau};
        solver.compute(parameters, boundaries_conditions, init_dist);
        for(const size_t step : std::ranges::iota_view{0u, 100001u}) {
            solver.calc_step(boundaries_conditions, right_part);
            if (step % 5 == 0) {
                auto solution = nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, parameters, solver.temperature()};
                solution.calc_flux();
                save_solution(solution, "./sol", step);
            }
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