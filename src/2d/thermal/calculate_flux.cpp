#include "influence_functions_2d.hpp"
#include "thermal/stationary_heat_equation_solver_2d.hpp"

namespace {

template<class T>
std::vector<T> read_solution(const std::string& path, const size_t nodes_count) {
    std::string str;
    std::ifstream file{path};
    std::vector<T> temperature(nodes_count);
    for(T& val : temperature) {
        file >> str;
        val = std::stod(str.substr(str.rfind(',') + 1));
    }
    return temperature;
}

}

int main(int argc, char** argv) {
    if(argc < 6) {
        std::cerr << "Input format [program name] <path to mesh> <path to solution> <r1> <r2> <p1>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(20);
        const double p1 = std::stod(argv[5]);
        const std::array<double, 2> r = {std::stod(argv[3]), std::stod(argv[4])};
        nonlocal::thermal::equation_parameters_2d<double, nonlocal::material_t::ORTHOTROPIC> eq_parameters;
        eq_parameters.lambda[0] = r[0] / std::max(r[0], r[1]);
        eq_parameters.lambda[1] = r[1] / std::max(r[0], r[1]);

        const auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int32_t>>(mesh);
        if (p1 < nonlocal::MAX_NONLOCAL_WEIGHT<double>) {
            mesh_proxy->find_neighbours(std::max(r[0], r[1]) + 0.015, nonlocal::mesh::balancing_t::MEMORY);
            double mean = 0;
            for(const size_t e : std::views::iota(size_t{0}, mesh->elements_count()))
                mean += mesh_proxy->neighbors(e).size();
            std::cout << "Average neighbours = " << mean / mesh->nodes_count() << std::endl;
        }

        const auto temperature = read_solution<double>(argv[2], mesh->nodes_count());
        nonlocal::thermal::solution<double, int32_t> T{mesh_proxy, eq_parameters,
                                                       p1, nonlocal::influence::polynomial_2d<double, 2, 1>{r},
                                                       temperature};
        T.calc_flux();
        std::cout << "Energy = " << T.calc_energy() << std::endl;
        T.save_as_vtk("heat.vtk");
        const auto& [X, Y] = T.flux();
        nonlocal::mesh::save_as_csv("T.csv", *mesh, T.temperature());
        nonlocal::mesh::save_as_csv("TX.csv", *mesh, X);
        nonlocal::mesh::save_as_csv("TY.csv", *mesh, Y);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}