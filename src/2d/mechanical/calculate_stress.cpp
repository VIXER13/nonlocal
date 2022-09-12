#include "influence_functions_2d.hpp"
#include "equilibrium_equation_2d.hpp"

namespace {

template<class T>
void save_raw_data(const std::string& path,
                   const nonlocal::mesh::mesh_2d<T>& msh,
                   const nonlocal::mechanical::mechanical_solution_2d<T, int>& sol) {
    std::ofstream
        u1{path + "/u1.csv"},
        u2{path + "/u2.csv"},
        eps11{path + "/eps11.csv"},
        eps22{path + "/eps22.csv"},
        eps12{path + "/eps12.csv"},
        sigma11{path + "/sigma11.csv"},
        sigma22{path + "/sigma22.csv"},
        sigma12{path + "/sigma12.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        u1      << msh.node(i)[0] << ',' << msh.node(i)[1] << "," << sol.displacement()[0][i] << std::endl;
        u2      << msh.node(i)[0] << ',' << msh.node(i)[1] << "," << sol.displacement()[0][i] << std::endl;
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[0][i] << std::endl;
        eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[1][i] << std::endl;
        eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[2][i] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[0][i] << std::endl;
        sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[1][i] << std::endl;
        sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[2][i] << std::endl;
    }
}

template<class T>
std::vector<T> read_solution(const std::string& path_u1, const std::string& path_u2, const size_t nodes_count) {
    std::ifstream u1{path_u1};
    std::ifstream u2{path_u2};
    std::vector<T> displacement(2 * nodes_count, 0);
    std::string str;
    for(const size_t i : std::views::iota(size_t{0}, nodes_count)) {
        u1 >> str;
        displacement[2 * i] = std::stod(str.substr(str.rfind(',') + 1));
        u2 >> str;
        displacement[2 * i + 1] = std::stod(str.substr(str.rfind(',') + 1));
    }
    return displacement;
}

}

int main(int argc, char** argv) {
    if(argc < 7) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1> <output_path> <path_to_u1> <path_to_u2>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);

        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            mesh_proxy->find_neighbours(r + 0.025, nonlocal::mesh::balancing_t::MEMORY);
        }

        nonlocal::mechanical::equation_parameters<double> parameters;
        parameters.nu = 0.3;
        parameters.E = 420;
        parameters.type = nonlocal::mechanical::calc_t::PLANE_STRESS;

        const std::vector<double> displacement = read_solution<double>(argv[5], argv[6], mesh->nodes_count());
        nonlocal::mechanical::mechanical_solution_2d<double, int64_t> sol{mesh_proxy, p1, bell, parameters, displacement};

        sol.calc_strain_and_stress();

        if (MPI_utils::MPI_rank() == 0) {
            std::cout << "Energy = " << sol.calc_energy() << std::endl;
            using namespace std::string_literals;
            sol.save_as_vtk(argv[4] + "/structural.vtk"s);
            save_raw_data(argv[4], *mesh, sol);
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}