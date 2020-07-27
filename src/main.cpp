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

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(5);
        omp_set_num_threads(4);

        static constexpr double r = 0.2, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};
        if(p1 < 1) {
            msh.find_elements_neighbors(r);
            size_t neighbors_count = 0;
            for(size_t i = 0; i < msh.elements_count(); ++i)
                neighbors_count += msh.element_neighbors(i).size();
            std::cout << "Average number of neighbors: " << double(neighbors_count) / msh.elements_count() << std::endl;
        }

        const nonlocal::structural::parameters<double> params = {.nu = 0.25, .E = 2.1e5};
        const auto u = nonlocal::structural::stationary(
            msh, 
            {
                {
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {
                    [](const std::array<double, 2>&) { return 0.1; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                },

                {
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                }
            },
            params, p1, bell
        );

        nonlocal::structural::save_as_vtk("test.vtk", msh, u);

        //msh.save_as_vtk("test.vtk", u);

        // msh.save_as_vtk(std::string{"test.vtk"},
        // {
        //     {std::string{"u_x"}, std::move(u_x)},
        //     {std::string{"u_y"}, std::move(u_y)}
        // });

        //raw_output("", msh, u);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}