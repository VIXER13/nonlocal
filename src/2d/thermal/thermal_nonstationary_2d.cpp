#include "influence_functions_2d.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"

namespace {

template<class T, class I>
void logger(const nonlocal::thermal::heat_equation_solution_2d<T, I>& solution, const std::string& path, const uintmax_t step) {
    std::cout << "step = " << step << std::endl;
    std::cout << "Energy = " << solution.calc_energy() << std::endl;
    solution.save_as_vtk(path + '/' + std::to_string(step) + ".vtk");
    const auto& [X, Y] = solution.flux();
    nonlocal::mesh::save_as_csv(path + '/' + std::to_string(step) + "T.csv",  solution.mesh_proxy()->mesh(), solution.temperature());
    if (!X.empty() && !Y.empty()) {
        nonlocal::mesh::save_as_csv(path + '/' + std::to_string(step) + "TX.csv", solution.mesh_proxy()->mesh(), X);
        nonlocal::mesh::save_as_csv(path + '/' + std::to_string(step) + "TY.csv", solution.mesh_proxy()->mesh(), Y);
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
        std::cout.precision(7);
        const double p1 = std::stod(argv[4]);
        const std::array<double, 2> r = {std::stod(argv[2]), std::stod(argv[3])};
        static const nonlocal::influence::polynomial_2d<double, 2, 1> influence_function(r);

        nonlocal::thermal::equation_parameters_2d<double, nonlocal::material_t::ORTHOTROPIC> eq_parameters;
        eq_parameters.lambda[0] = r[0] / std::max(r[0], r[1]);
        eq_parameters.lambda[1] = r[1] / std::max(r[0], r[1]);
        eq_parameters.rho = 1;
        eq_parameters.c = 1;
        const double tau = 0.01;

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            const double time = omp_get_wtime();
            mesh_proxy->find_neighbours(std::max(r[0], r[1]) + 0.015, nonlocal::mesh::balancing_t::MEMORY);
            std::cout << "find neighbours: " << omp_get_wtime() - time << std::endl;
            double mean = 0;
            for(const size_t e : std::views::iota(size_t{0}, mesh->elements_count()))
                mean += mesh_proxy->neighbors(e).size();
            std::cout << "Average neighbours = " << mean / mesh->elements_count() << std::endl;
        }

        const std::unordered_map<std::string, nonlocal::nonstationary_boundary_2d_t<nonlocal::thermal::boundary_condition_t, double, 1>>
            boundary_conditions = {
                {   "Right",
                    {   nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return 1; },
                    }
                },
                {   "Left",
                    {   nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return -1; },
                    }
                }
        };

        nonlocal::thermal::nonstationary_heat_equation_solver_2d<double, int32_t, int32_t> solver{mesh_proxy, tau};
        solver.compute(eq_parameters, nonlocal::boundary_type(boundary_conditions),
                       [](const std::array<double, 2>&) constexpr noexcept { return 0; },
                       p1, influence_function);
        logger(nonlocal::thermal::heat_equation_solution_2d<double, int32_t>{mesh_proxy}, argv[5], 0);
        for(const uintmax_t step : std::views::iota(1, 101)) {
            solver.calc_step(eq_parameters.alpha, boundary_conditions,
                [](const double t, const std::array<double, 2>& x) constexpr noexcept { return 0; });
            auto solution = nonlocal::thermal::heat_equation_solution_2d<double, int32_t>{mesh_proxy, eq_parameters, p1, influence_function, solver.temperature()};
            solution.calc_flux();
            logger(solution, argv[5], step);
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