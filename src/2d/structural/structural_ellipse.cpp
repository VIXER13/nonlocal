#include "solver_2d/influence_functions_2d.hpp"
#include "solver_2d/structural/structural_solver.hpp"

namespace {

template<class T>
void save_raw_data(const std::string& path,
                   const nonlocal::mesh::mesh_2d<T>& msh,
                   const nonlocal::structural::solution<T, int>& sol) {
    std::ofstream eps11{path + "/eps11.csv"},
            eps22{path + "/eps22.csv"},
            eps12{path + "/eps12.csv"},
            sigma11{path + "/sigma11.csv"},
            sigma22{path + "/sigma22.csv"},
            sigma12{path + "/sigma12.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[0][i] << std::endl;
        eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[1][i] << std::endl;
        eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[2][i] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [0][i] << std::endl;
        sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [1][i] << std::endl;
        sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [2][i] << std::endl;
    }
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
        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell{r};

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999)
            mesh_proxy->find_neighbours(r + 0.05, nonlocal::mesh::balancing_t::MEMORY);

        nonlocal::structural::structural_solver<double, int, int> fem_sol{mesh_proxy};
        nonlocal::structural::equation_parameters<double> parameters;
        parameters.nu = 0.3;
        parameters.E = 21;
        parameters.p1 = p1;

        auto sol = fem_sol.stationary(parameters,
            { // Граничные условия
                {   // Left
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return -1; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Right
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 1; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Hor
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Ver
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                }
            },
            {}, // Правая часть
            bell // Функция влияния
        );

        sol.calc_strain_and_stress();

        if (MPI_utils::MPI_rank() == 0) {
            std::cout << "Energy = " << sol.calc_energy() << std::endl;
            using namespace std::string_literals;
            sol.save_as_vtk(argv[4] + "/structural.vtk"s);
            save_raw_data(argv[4], *mesh, sol);
        }
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