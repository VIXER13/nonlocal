#include "solver_2d/influence_functions_2d.hpp"
#include "solver_2d/thermal/heat_equation_solver_2d.hpp"

int main(int argc, char** argv) {
#ifdef MPI_USE
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
#endif

    if(argc < 5) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1> <save_path>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);
        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);
        nonlocal::heat::solver_parameters<double> sol_parameters;
        sol_parameters.save_path = argv[4];
        sol_parameters.time_interval = {0, 5};
        sol_parameters.steps = 10000;
        sol_parameters.save_freq = 50;
        sol_parameters.save_vtk = false;
        sol_parameters.save_csv = true;
        sol_parameters.calc_energy = true;

        nonlocal::heat::equation_parameters<double, nonlocal::heat::material_t::ISOTROPIC> eq_parameters;
        eq_parameters.lambda[0] = 1;
        eq_parameters.rho = 1;
        eq_parameters.c = 1;
        eq_parameters.p1 = p1;

        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            mesh_proxy->find_neighbours(r + 0.05, mesh::balancing_t::MEMORY); // 0.05 это некая гарантированная добавка
            // Нужна, чтобы все квадратурные узлы, которые попадают под зону влияния были учтены.
            double mean = 0;
            for(size_t e = 0; e < mesh->elements_count(); ++e)
                mean += mesh_proxy->neighbors(e).size();
            std::cout << "Average neighbours = " << mean << std::endl;
        }
        nonlocal::heat::heat_equation_solver_2d<double, int, int> fem_sol{mesh_proxy};

        fem_sol.nonstationary(
            sol_parameters, eq_parameters,
            { // Граничные условия
                {   // Down
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return -1; },
                },

                {   // Right
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                {   // Up
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 1; },
                },

                {   // Left
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                }
            },
            [](const std::array<double, 2>&) { return 0; }, // Начальные условия
            [](const std::array<double, 2>&) { return 0; }, // Правая часть
            bell // Функция влияния
        );
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