#include "solvers/influence_functions.hpp"
#include "solvers/heat_equation_solver.hpp"

namespace {

void save_raw_data(const mesh::mesh_2d<double>& msh, const Eigen::Matrix<double, Eigen::Dynamic, 1>& T) {
    std::ofstream Tout{"T.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i)
        Tout << msh.node(i)[0] << ',' << msh.node(i)[1] << ',' << T[i] << '\n';
}

}

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);
        omp_set_num_threads(4);

        static constexpr double r = 0.2, p1 = 1;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};

        const auto& e = msh.element_2d(0);
        
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            std::cout << e->weight(q) << ' ';
        std::cout << std::endl;
        std::cout << std::endl;
        
        for(size_t i = 0; i < e->nodes_count(); ++i) {
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                std::cout << e->qN(i, q) << ' ';
            std::cout << std::endl;
        }
        std::cout << std::endl;

        for(size_t i = 0; i < e->nodes_count(); ++i) {
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                std::cout << e->qNxi(i, q) << ' ';
            std::cout << std::endl;
        }

        std::cout << std::endl;

        for(size_t i = 0; i < e->nodes_count(); ++i) {
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                std::cout << e->qNeta(i, q) << ' ';
            std::cout << std::endl;
        }

        std::cout << std::endl;



        nonlocal::heat::heat_equation_solver<double, int> fem_sol{msh};

        const auto T = fem_sol.stationary(
            { // Граничные условия
                {   // Down
                        [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                        nonlocal::heat::boundary_t::FLOW
                },

                {   // Right
                        [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                        nonlocal::heat::boundary_t::FLOW
                },

                {   // Up
                        [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                        nonlocal::heat::boundary_t::FLOW
                },

                {   // Left
                        [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                        nonlocal::heat::boundary_t::FLOW
                }
            }, 
            [](const std::array<double, 2>&) { return 0; }, // Правая часть
            r,  // Радиус влияния 
            p1, // Вес
            bell // Функция влияния
        );

        fem_sol.save_as_vtk("heat.vtk", T);
        save_raw_data(msh, T);

    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}