#include "influence_functions_2d.hpp"
#include "equilibrium_equation_2d.hpp"

namespace {

template<class T>
void save_raw_data(const std::string& path,
                   const nonlocal::mesh::mesh_2d<T>& msh,
                   const nonlocal::mechanical::mechanical_solution_2d<T, int>& sol) {
    std::ofstream eps11{path + "/eps11.csv"},
            eps22{path + "/eps22.csv"},
            eps12{path + "/eps12.csv"},
            sigma11{path + "/sigma11.csv"},
            sigma22{path + "/sigma22.csv"},
            sigma12{path + "/sigma12.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[0][i] << std::endl;
        eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[1][i] << std::endl;
        eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strain()[2][i] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[0][i] << std::endl;
        sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[1][i] << std::endl;
        sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress()[2][i] << std::endl;
    }
}

double f1(const std::array<double, 2>& x) { return 1; }
double f2(const std::array<double, 2>& x) { return x[1] < 0.5 ? x[1] : 1 - x[1]; }
double f3(const std::array<double, 2>& x) { return f2(x) - 0.25; }
double f4(const std::array<double, 2>& x) { return 0.5 - f2(x); }
double f5(const std::array<double, 2>& x) { return 0.25 - f2(x); }
double f6(const std::array<double, 2>& x) { return x[1] < 0.5 ? 0.5 - x[1] : x[1] - 0.5; }
double f7(const std::array<double, 2>& x) { return x[1] < 0.45 || x[1] > 0.55 ? 0 : 10; }

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
        std::cout.precision(7);

        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999) {
            mesh_proxy->find_neighbours(r + 0.025, nonlocal::mesh::balancing_t::MEMORY);
        }

        nonlocal::mechanical::equation_parameters<double> parameters;
        parameters.poissons_ratio = 0.3;
        parameters.young_modulus = 420;
        parameters.task = nonlocal::mechanical::task_2d_t::PLANE_STRESS;

//        parameters.thermoelasticity = true;
//        parameters.delta_temperature.resize(mesh->nodes_count());
//        parameters.alpha = 1;
//        for(size_t i = 0; i < mesh->nodes_count(); ++i) {
//            const auto& node = mesh->node(i);
//            parameters.delta_temperature[i] = node[0] + node[1];
//        }

        auto sol = nonlocal::mechanical::equilibrium_equation<double, int, int64_t>(parameters, mesh_proxy,
            { // Граничные условия
                {   "Right",
                    {
                    nonlocal::mechanical::boundary_condition_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 1; },
                    nonlocal::mechanical::boundary_condition_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                    }
                },

//                {   "Left",
//                    {
//                    nonlocal::mechanical::boundary_condition_t::PRESSURE,
//                    [](const std::array<double, 2>& x) { return -1; },
//                    nonlocal::mechanical::boundary_condition_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 0; }
//                    }
//                },

                {   "Horizontal",
                    {
                        nonlocal::mechanical::boundary_condition_t::PRESSURE,
                        [](const std::array<double, 2>&) { return 0; },
                        nonlocal::mechanical::boundary_condition_t::DISPLACEMENT,
                        [](const std::array<double, 2>&) { return 0; }
                    }
                },

                {   "Vertical",
                    {
                        nonlocal::mechanical::boundary_condition_t::DISPLACEMENT,
                        [](const std::array<double, 2>&) { return 0; },
                        nonlocal::mechanical::boundary_condition_t::PRESSURE,
                        [](const std::array<double, 2>&) { return 0; }
                    }
                }
            },
            [](const std::array<double, 2>&) { return std::array<double, 2>{}; }, // Правая часть
            p1,
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