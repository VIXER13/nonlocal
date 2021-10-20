#include "influence_functions_2d.hpp"
#include "structural_solver.hpp"

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

template<class T>
std::vector<T> read_displacement(const std::string& path_to_u1, const std::string& path_to_u2, const size_t size) {
    T val = 0;
    std::vector<T> u;
    std::ifstream file_u1{path_to_u1}, file_u2{path_to_u2};

    std::string str1, str2;
    while (std::getline(file_u1, str1) && std::getline(file_u2, str2)) {
        for(char& c : str1)
            if (c == ',')
                c = ' ';
        std::stringstream stream1{str1};
        stream1 >> val >> val >> val;
        u.push_back(val);

        std::stringstream stream2{str2};
        for(char& c : str2)
            if (c == ',')
                c = ' ';
        stream2 >> val >> val >> val;
        u.push_back(val);
    }
    return u;
}

}

int main(int argc, char** argv) {
#if MPI_USE
    MPI_Init(&argc, &argv);
#endif

    if(argc < 7) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1> <output_path> <u1_path> <u2_path>";
#if MPI_USE
        MPI_Finalize();
#endif
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

        std::cout << argv[5] << std::endl;
        std::cout << argv[6] << std::endl;
        const auto displacement = read_displacement<double>(argv[5], argv[6], mesh->nodes_count());
        std::cout << mesh->nodes_count() << std::endl;
        std::cout << displacement.size() << std::endl;

        nonlocal::structural::equation_parameters<double> parameters;
        parameters.nu = 0.3;
        parameters.E = 420;
        parameters.p1 = p1;
        parameters.type = nonlocal::structural::calc_t::PLANE_STRESS;
        nonlocal::structural::solution<double, int> sol{mesh_proxy, parameters, bell, displacement};

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