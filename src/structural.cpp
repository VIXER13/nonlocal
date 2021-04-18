#include "influence_functions.hpp"
#include "structural_solver.hpp"

namespace {

template<class T>
void save_raw_data(const mesh::mesh_2d<T>& msh,
                   const nonlocal::structural::solution<T, int>& sol) {
    std::ofstream eps11{"eps11.csv"},
                  eps22{"eps22.csv"},
                  eps12{"eps12.csv"},
                  sigma11{"sigma11.csv"},
                  sigma22{"sigma22.csv"},
                  sigma12{"sigma12.csv"};
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
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 5) {
        std::cerr << "Input format [program name] <path to mesh> <num_threads> <r> <p1>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);

        omp_set_num_threads(std::stoi(argv[2]));
        //mkl_set_num_threads(std::stoi(argv[2]));

        const double r = std::stod(argv[3]), p1 = std::stod(argv[4]);
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);

        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999)
            mesh_proxy->find_neighbours(r, mesh::balancing_t::MEMORY);

        nonlocal::structural::structural_solver<double, int, long long> fem_sol{mesh_proxy};
        nonlocal::structural::calculation_parameters<double> parameters;
        parameters.nu = 0.3;
        parameters.E = 21;

        auto sol = fem_sol.stationary(parameters,
            { // Граничные условия
                {   // Down
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Right
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 1; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Up
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; }
                },

                {   // Left
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>& x) { return 0; }
                }
            },
            {}, // Правая часть
            p1, // Вес
            bell // Функция влияния
        );

        sol.calc_strain_and_stress();

        if (mesh_proxy->rank() == 0) {
            std::cout << "Energy = " << sol.calc_energy() << std::endl;
            sol.save_as_vtk("structural.vtk");
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        PetscFinalize();
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        PetscFinalize();
        return EXIT_FAILURE;
    }

    PetscFinalize();
    return EXIT_SUCCESS;
}