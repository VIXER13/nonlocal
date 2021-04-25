#include <petsc.h>
#include "influence_functions.hpp"
#include "heat_equation_solver.hpp"

namespace {

    void save_raw_data(const std::shared_ptr<mesh::mesh_2d<double>>& mesh,
                       const std::vector<double>& T,
                       const std::array<std::vector<double>, 2>& gradient) {
        std::ofstream Tout{"T.csv"};
        std::ofstream Tx{"Tx.csv"};
        std::ofstream Ty{"Ty.csv"};
        for(size_t i = 0; i < mesh->nodes_count(); ++i) {
            Tout << mesh->node(i)[0] << ',' << mesh->node(i)[1] << ',' << T[i] << '\n';
            Tx   << mesh->node(i)[0] << ',' << mesh->node(i)[1] << ',' << gradient[0][i] << '\n';
            Ty   << mesh->node(i)[0] << ',' << mesh->node(i)[1] << ',' << gradient[1][i] << '\n';
        }
    }

}

int main(int argc, char** argv) {
#ifdef MPI_USE
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
#endif

    if(argc < 4) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);
        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        nonlocal::heat::calculation_parameters<double, nonlocal::heat::calc_type::ISOTROPIC> parameters;
        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<mesh::mesh_proxy<double, int>>(mesh);
        nonlocal::heat::heat_equation_solver<double, int, int> fem_sol{mesh_proxy};
        if (p1 < 0.999)
            mesh_proxy->find_neighbours(r, mesh::balancing_t::MEMORY);

        fem_sol.nonstationary(
            ".",
            parameters, 0.0001, 1000,
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
            p1, // Вес
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