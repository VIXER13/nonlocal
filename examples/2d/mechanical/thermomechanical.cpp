#include "influence_functions_2d.hpp"
#include "equilibrium_equation_2d.hpp"
#include "stationary_heat_equation_solver_2d.hpp"

namespace nonlocal {

template<class T, class I, class Influence>
nonlocal::thermal::heat_equation_solution_2d<T, I> thermal_equation(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                    const T local_weight, const Influence& influence) {
    static constexpr T conductivity = 1;
    static constexpr T capacity = 1;
    static constexpr T density = 1;
    const nonlocal::thermal::parameters_2d<T> parameters = {{
        "DEFAULT",
        {
            .model = {
                .influence = influence,
                .local_weight = local_weight
            },
            .physical = std::make_shared<thermal::parameter_2d<T, coefficients_t::CONSTANTS>>(
                conductivity, capacity, density
            )
        }
    }};
    nonlocal::thermal::thermal_boundaries_conditions_2d<T> conditions;
    conditions["Left"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(-1);
    conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(1);
    auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I, T, I>(
        mesh, parameters, conditions,
        [](const std::array<T, 2>&) constexpr noexcept { return 0; }
    );
    solution.calc_flux();
    return solution;
}

template<class T, class I, class Influence>
nonlocal::mechanical::mechanical_solution_2d<T, I> equilibrium_equation(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                        const std::vector<T>& temperature,
                                                                        const T local_weight, const Influence& influence) {
    static constexpr T E = 400;
    static constexpr T nu = 0.3;
    static constexpr T alpha = 2.5e-3;
    const nonlocal::mechanical::mechanical_parameters_2d<T> parameters = {
        .materials = {{
            "DEFAULT",
            {
                .model = {
                    .influence = influence,
                    .local_weight = local_weight
                },
                .physical = {
                    .young_modulus = E,
                    .poissons_ratio = nu,
                    .thermal_expansion = alpha
                }
            }
        }},
        .delta_temperature = temperature
    };
    nonlocal::mechanical::mechanical_boundaries_conditions_2d<T> conditions;
    conditions["Left"] = {
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0),
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0)
    };
    conditions["Right"] = {
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0),
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0)
    };
    conditions["Horizontal"] = {
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0),
        std::make_unique<nonlocal::mechanical::displacement_2d<T>>(0)
    };
    conditions["Vertical"] = {
        std::make_unique<nonlocal::mechanical::displacement_2d<T>>(0),
        std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0)
    };
    auto solution = nonlocal::mechanical::equilibrium_equation<I>(
        mesh, parameters, conditions,
        [](const std::array<T, 2>&) constexpr noexcept { return std::array<T, 2>{}; }
    );
    solution.calc_strain_and_stress();
    return solution;
}

}

int main(int argc, char** argv) {
#if MPI_USE
    MPI_Init(&argc, &argv);
#endif

    if(argc < 5) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1> <output_path>";
#if MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    int exit_code = EXIT_SUCCESS;
    try {
        using T = double;
        using I = int64_t;

        const std::unordered_map<std::string, T> r = {{"DEFAULT", std::stod(argv[2])}};
        const double p1 = std::stod(argv[3]);
        const nonlocal::influence::polynomial_2d<double, 2, 1> influence{r.at("DEFAULT")};
        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
        if (nonlocal::theory_type(p1) == nonlocal::theory_t::NONLOCAL)
            mesh->find_neighbours(r);
        
        const auto thermal_sol = nonlocal::thermal_equation(mesh, p1, influence);
        nonlocal::mesh::utils::save_as_csv("temperature.csv", mesh->container(), thermal_sol.temperature());
        nonlocal::mesh::utils::save_as_csv("flux_x.csv", mesh->container(), thermal_sol.flux()[0]);
        nonlocal::mesh::utils::save_as_csv("flux_y.csv", mesh->container(), thermal_sol.flux()[1]);
        thermal_sol.save_as_vtk("thermal.vtk");
        
        const auto mech_sol = nonlocal::equilibrium_equation(mesh, thermal_sol.temperature(), p1, influence);
        nonlocal::mesh::utils::save_as_csv("displacement_x.csv", mesh->container(), mech_sol.displacement()[0]);
        nonlocal::mesh::utils::save_as_csv("displacement_y.csv", mesh->container(), mech_sol.displacement()[1]);
        nonlocal::mesh::utils::save_as_csv("strain_11.csv", mesh->container(), mech_sol.strain()[0]);
        nonlocal::mesh::utils::save_as_csv("strain_22.csv", mesh->container(), mech_sol.strain()[1]);
        nonlocal::mesh::utils::save_as_csv("strain_12.csv", mesh->container(), mech_sol.strain()[2]);
        nonlocal::mesh::utils::save_as_csv("stress_11.csv", mesh->container(), mech_sol.stress()[0]);
        nonlocal::mesh::utils::save_as_csv("stress_22.csv", mesh->container(), mech_sol.stress()[1]);
        nonlocal::mesh::utils::save_as_csv("stress_12.csv", mesh->container(), mech_sol.stress()[2]);
        mech_sol.save_as_vtk("mechanical.vtk");
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
#if MPI_USE
        MPI_Finalize();
#endif
        exit_code = EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
#if MPI_USE
        MPI_Finalize();
#endif
        exit_code = EXIT_FAILURE;
    }

#if MPI_USE
    MPI_Finalize();
#endif

    return exit_code;
}