#include "solvers/influence_functions.hpp"
#include "solvers/heat_equation_solver.hpp"

template<class Type, class Index, class Vector>
static void raw_output(const std::string& path, const mesh::mesh_2d<Type, Index>& mesh, const Vector& T)
{
    std::ofstream fout(path + std::string{"T.csv"});
    fout.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << mesh.node(i)[0] << "," << mesh.node(i)[1] << "," << T[i] << std::endl;
}

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(5);
        omp_set_num_threads(4);

        static constexpr double r = 0.2, p1 = 1;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
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
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },
                
                {
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },

                {
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },

                {
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                }
            },
            [](const std::array<double, 2>&){ return -4; },
            p1, bell, 1.
        );

        std::cout << "Energy: " << nonlocal::heat::integrate_solution(msh, T) << std::endl;

        //msh.save_as_vtk("test.vtk", T);
        raw_output("test.csv", msh, T);

        nonlocal::heat::nonstationary("results/nonstationary/",
            msh, 0.0001, 1000,
            {
                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                },
                
                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                },

                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                },

                {
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::FLOW
                }
            },
            [](const std::array<double, 2>& x) { return 1.2 * (1. - metamath::power<2>(x[0] - 0.5) -
                                                                    metamath::power<2>(x[1] - 0.5)); },
            [](const std::array<double, 2>&) { return 0; },
            p1, bell, 1
        );
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}