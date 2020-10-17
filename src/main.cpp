#include "solvers/influence_functions.hpp"
#include "solvers/structural_solver.hpp"

template<class Type, class Index, class Vector>
static void raw_output(const std::string& path, const mesh::mesh_2d<Type, Index>& mesh, const Vector& T) {
    std::ofstream fout(path + std::string{"T.csv"});
    fout.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << mesh.node(i)[0] << "," << mesh.node(i)[1] << "," << T[i] << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "test" << std::endl;
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    std::cout << "test" << std::endl;

    try {
        std::cout.precision(16);
        omp_set_num_threads(4);

        static constexpr double r = 0.05, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};

        std::cout << "test" << std::endl;
        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 42};
        nonlocal::structural::structural_solver<double, int> fem_sol{msh, params};

        const auto U = fem_sol.stationary(
            {
                {},

                {
                    [](const std::array<double, 2>&) { return 1; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {},

                {
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                },

                {}, {}
            },
            r, p1, bell
        );

        const auto [strain, stress] = fem_sol.strains_and_stress(U, p1, bell);
        std::cout << "Energy U = " << fem_sol.calc_energy(strain, stress) << std::endl;
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}