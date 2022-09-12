#include "influence_functions_2d.hpp"
#include "equilibrium_equation_2d.hpp"

namespace {
template<class T, class I, class Influence_Function>
nonlocal::mechanical::mechanical_solution_2d<T, I> mechanic(
    const std::shared_ptr<nonlocal::mesh::mesh_proxy<T, I>>& mesh_proxy,
    const std::vector<T>& delta_temperature,
    const T p1,
    const Influence_Function& influence_function) {
    nonlocal::mechanical::equation_parameters<double> params;
    params.task = nonlocal::mechanical::task_2d_t::PLANE_STRESS;
    params.young_modulus = 210e9;
    params.poissons_ratio = 0.3;
    params.thermal_expansion = 11.9e-6;
    params.is_thermoelasticity = true;
    params.delta_temperature = delta_temperature;

    const std::unordered_map<std::string, nonlocal::stationary_boundary_2d_t<nonlocal::mechanical::boundary_condition_t, T, 2>>
    boundary_condition = { // Граничные условия
        {   "Right",
            {
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>& x) { return 0; },
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; }
            }
        },

        {  "Left",
            {
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>& x) { return 0; },
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; }
            }
        },

        {   "Down",
            {
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>& x) { return 0; },
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; }
            }
        },

        {  "Up",
            {
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>& x) { return 0; },
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; }
            }
        },

        {   "Horizontal",
            {
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; },
                nonlocal::mechanical::boundary_condition_t::DISPLACEMENT,
                [](const std::array<double, 2>&) { return 0; }
            }
        },

        {   "Vertical",
            {
                nonlocal::mechanical::boundary_condition_t::DISPLACEMENT,
                [](const std::array<double, 2>&) { return 0; },
                nonlocal::mechanical::boundary_condition_t::PRESSURE,
                [](const std::array<double, 2>&) { return 0; }
            }
        }
    };

    return nonlocal::mechanical::equilibrium_equation<T, I, I>(
        params, mesh_proxy, boundary_condition,
        [](const std::array<T, 2> &x) { return std::array<double, 2>{}; }, // Правая часть
        p1, influence_function
    );
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

    try {
        std::cout.precision(7);

        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            const double time = omp_get_wtime();
            mesh_proxy->find_neighbours(r + 0.015, nonlocal::mesh::balancing_t::MEMORY);
            std::cout << "find neighbours time: " << omp_get_wtime() - time << std::endl;
            double mean = 0;
            for(const size_t e : std::views::iota(size_t{0}, mesh->elements_count()))
                mean += mesh_proxy->neighbors(e).size();
            std::cout << "Average neighbours = " << mean / mesh->nodes_count() << std::endl;
        }

        const std::vector<double> temperature(mesh->nodes_count(), 1);
        auto mechanical_sol = mechanic(mesh_proxy, temperature, p1, bell);
        mechanical_sol.calc_strain_and_stress();

        using namespace std::literals;
        mechanical_sol.save_as_vtk(argv[4] + "/mechanical.vtk"s);
        nonlocal::mesh::save_as_csv(argv[4] + "/u1.csv"s, *mesh, mechanical_sol.displacement()[0]);
        nonlocal::mesh::save_as_csv(argv[4] + "/u2.csv"s, *mesh, mechanical_sol.displacement()[1]);
        nonlocal::mesh::save_as_csv(argv[4] + "/eps11.csv"s, *mesh, mechanical_sol.strain()[0]);
        nonlocal::mesh::save_as_csv(argv[4] + "/eps22.csv"s, *mesh, mechanical_sol.strain()[1]);
        nonlocal::mesh::save_as_csv(argv[4] + "/eps12.csv"s, *mesh, mechanical_sol.strain()[2]);
        nonlocal::mesh::save_as_csv(argv[4] + "/sigma11.csv"s, *mesh, mechanical_sol.stress()[0]);
        nonlocal::mesh::save_as_csv(argv[4] + "/sigma22.csv"s, *mesh, mechanical_sol.stress()[1]);
        nonlocal::mesh::save_as_csv(argv[4] + "/sigma12.csv"s, *mesh, mechanical_sol.stress()[2]);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
#if MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
#if MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

#if MPI_USE
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}