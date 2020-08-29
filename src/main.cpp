#include "solvers/influence_functions.hpp"
#include "solvers/static_equation_solver.hpp"

template<class Type, class Index, class Vector>
static void raw_output(const std::string& path, const mesh::mesh_2d<Type, Index>& mesh, const Vector& u) {
    std::ofstream fu1{path + std::string{"u1.csv"}},
                  fu2{path + std::string{"u2.csv"}};
    fu1.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fu1 << mesh.node(i)[0] << "," << mesh.node(i)[1] << "," << u[2*i] << std::endl;

    fu2.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fu2 << mesh.node(i)[0] << "," << mesh.node(i)[1] << "," << u[2*i+1] << std::endl;
}

template<class Type>
static std::tuple<std::string, std::string, Type, Type>
    parse_parameters(const std::string& path) {
    std::ifstream parameters{path};
    std::string in, out;
    Type r = 0, p1 = 0;
    parameters >> in >> out >> r >> p1;
    return {std::move(in), std::move(out), r, p1};
}

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to init>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(5);
        omp_set_num_threads(4);

        const auto [in_file, out_file, r, p1] = parse_parameters<double>(argv[1]);
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> mesh{in_file};

        if(p1 < 1) {
            mesh.find_elements_neighbors(1.05*r);
            size_t neighbors_count = 0;
            for(size_t i = 0; i < mesh.elements_count(); ++i)
                neighbors_count += mesh.element_neighbors(i).size();
            std::cout << "Average number of neighbors: " << double(neighbors_count) / mesh.elements_count() << std::endl;
        }

        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 2.1e5};
        const auto u = nonlocal::structural::stationary(
            mesh, 
            {
                {   // Down
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return -10000; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {}, // Right

                {   // Horizontal
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                },

                {}, // Left

                {   // Vertical
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {}, // EllipseUpDown

                {   // Up
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 10000; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {}  // EllipseLeftRight
            },
            params, p1, bell
        );

        if(p1 < 1)
            mesh.find_nodes_neighbors(1.05*r);
        auto [strain, stress] = nonlocal::structural::strains_and_stress(mesh, params, u, p1, bell);

        std::cout.precision(10);
        std::cout << "Energy: " <<  nonlocal::structural::calc_energy(mesh, strain, stress) << std::endl;

        nonlocal::structural::save_as_vtk(out_file, mesh, u, strain, stress);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}