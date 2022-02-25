#include "influence_functions_2d.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"

namespace {

template<class Vector>
void logger(const std::string& path,
            const Vector& temperature,
            const nonlocal::thermal::equation_parameters_2d<double, nonlocal::material_t::ORTHOTROPIC>& eq_parameters,
            const nonlocal::mesh::mesh_proxy<double, int>& mesh_proxy, const uintmax_t step) {
    std::cout << "step = " << step << std::endl;

    auto flux = mesh_proxy.gradient(temperature);
    for(size_t comp = 0; comp < flux.size(); ++comp)
        for(double& val : flux[comp])
            val *= eq_parameters.lambda[comp];

    nonlocal::mesh::save_as_vtk(path + "/"  + std::to_string(step) + ".vtk", mesh_proxy.mesh(), temperature);
    nonlocal::mesh::save_as_vtk(path + "/X" + std::to_string(step) + ".vtk", mesh_proxy.mesh(), flux[0]);
    nonlocal::mesh::save_as_vtk(path + "/Y" + std::to_string(step) + ".vtk", mesh_proxy.mesh(), flux[1]);

    nonlocal::mesh::save_as_csv(path + "/"  + std::to_string(step) + ".csv", mesh_proxy.mesh(), temperature);
    nonlocal::mesh::save_as_csv(path + "/X" + std::to_string(step) + ".csv", mesh_proxy.mesh(), flux[0]);
    nonlocal::mesh::save_as_csv(path + "/Y" + std::to_string(step) + ".csv", mesh_proxy.mesh(), flux[1]);

    std::cout << "Energy = " << mesh_proxy.integrate_solution(temperature) << std::endl;
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
        std::cout << "r[0] = " << r[0] << std::endl;
        std::cout << "r[1] = " << r[1] << std::endl;
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);

        nonlocal::thermal::equation_parameters_2d<double, nonlocal::material_t::ORTHOTROPIC> eq_parameters;
        eq_parameters.lambda[0] = r[0] / std::max(r[0], r[1]);
        eq_parameters.lambda[1] = r[1] / std::max(r[0], r[1]);
        eq_parameters.rho = 1;
        eq_parameters.c = 1;
        const double tau = 0.01;

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            mesh_proxy->find_neighbours(std::max(r[0], r[1]) + 0.05, nonlocal::mesh::balancing_t::MEMORY); // 0.05 это некая гарантированная добавка
            // Нужна, чтобы все квадратурные узлы, которые попадают под зону влияния были учтены.
            double mean = 0;
            for(const size_t e : std::views::iota(size_t{0}, mesh->elements_count()))
                mean += mesh_proxy->neighbors(e).size();
            std::cout << "Average neighbours = " << mean / mesh->elements_count() << std::endl;
        }

        const std::unordered_map<std::string, nonlocal::nonstationary_boundary_2d_t<nonlocal::thermal::boundary_condition_t, double, 1>>
            boundary_conditions = {
                {   "Down",
                    {
                        nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return 0; },
                    }
                },
                {   "Right",
                    {
                        nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return 1; },
                    }
                },
                {   "Up",
                    {
                        nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return 0; },
                    }
                },
                {   "Left",
                    {
                        nonlocal::thermal::boundary_condition_t::FLUX,
                        [](const double t, const std::array<double, 2>& x) { return -1; },
                    }
                }
        };

        nonlocal::thermal::nonstationary_heat_equation_solver_2d<double, int32_t, int32_t> solver{mesh_proxy, tau};
        solver.compute(eq_parameters, nonlocal::boundary_type(boundary_conditions),
                       [](const std::array<double, 2>&) constexpr noexcept { return 0; },
                       p1, bell);
        logger(argv[5], solver.temperature(), eq_parameters, *mesh_proxy, 0);
        for(const uintmax_t step : std::views::iota(1, 101)) {
            solver.calc_step(boundary_conditions, [](const double t, const std::array<double, 2>& x) constexpr noexcept { return 0; });
            logger(argv[5], solver.temperature(), eq_parameters, *mesh_proxy, step);
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