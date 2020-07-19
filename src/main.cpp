#include "solvers/influence_functions.hpp"
#include "solvers/heat_equation_solver.hpp"

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        static constexpr double r = 0.2, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 2> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};
        if(p1 < 1) {
            msh.find_elements_neighbors(r);
            size_t neighbors_count = 0;
            for(size_t i = 0; i < msh.elements_count(); ++i)
                neighbors_count += msh.element_neighbors(i).size();
            std::cout << "Average number of neighbors: " << double(neighbors_count) / msh.elements_count() << std::endl;
        }

        const auto T = nonlocal::heat::stationary(msh,
            {
                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },
                
                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                },

                {
                    [](const std::array<double, 2>&) { return 1; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },

                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                }
            },
            [](const std::array<double, 2>&){ return 0; },
            p1, bell
        );

        msh.save_as_vtk("test.vtk", T);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}