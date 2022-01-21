#include "influence_functions_2d.hpp"
#include "structural_solver.hpp"
#include "heat_equation_solver_2d.hpp"

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
        std::cout.precision(15);

        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double>>(argv[1]);

        const double r = std::stod(argv[2]);//, p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial_2d<double, 2, 1> bell(r);

        //for(const double p1 : {1., 2./3., 0.5, 1./3.}) {
        for(const double p1 : {2./3.}) {
            std::cout << "p1 = " << p1 << std::endl << std::endl;
            //for(const double s : {-1/3., -2/9., -1/9., 0., 1/9., 2/9., 1/3., 4/9., 5/9., 2/3.}) {
            for(const double s : {2/9.}) {
                auto *const el = dynamic_cast<metamath::finite_element::element_2d_integrate<double, metamath::finite_element::quadratic_serendipity>*>(mesh->element_2d(0).get());
                std::cout << "s = " << s << std::endl;
                el->set_parameter(s);
                el->set_quadrature(metamath::finite_element::quadrature_1d<double, metamath::finite_element::gauss3>{},
                                   metamath::finite_element::quadrature_1d<double, metamath::finite_element::gauss3>{});

                auto mesh_proxy = std::make_shared<nonlocal::mesh::mesh_proxy<double, int>>(mesh);
                if (p1 < 0.999) {
                    mesh_proxy->find_neighbours(r + 0.015, nonlocal::mesh::balancing_t::MEMORY);
                }

                nonlocal::heat::heat_equation_solver_2d<double, int, int> fem_sol{mesh_proxy};


            nonlocal::heat::equation_parameters<double, nonlocal::heat::material_t::ISOTROPIC> eq_parameters;
        //eq_parameters.lambda[0] = r[0] / std::max(r[0], r[1]);
        //eq_parameters.lambda[1] = r[1] / std::max(r[0], r[1]);
        eq_parameters.lambda = {1};
        eq_parameters.p1 = p1;
        //eq_parameters.r  = r;

        auto T = fem_sol.stationary(eq_parameters,
            { // Граничные условия
                {   "Right",
                    {
                        nonlocal::heat::boundary_t::FLOW,
                        [](const std::array<double, 2>& x) { return 1; },
                    }
                },

                {   "Horizontal",
                    {
                        nonlocal::heat::boundary_t::FLOW,
                        [](const std::array<double, 2>& x) { return 0; },
                    }
                },

                {
                    "Left",
                    {
                        nonlocal::heat::boundary_t::FLOW,
                        [](const std::array<double, 2>& x) { return -1; },
                    }
                },

                {   "Vertical",
                    {
                        nonlocal::heat::boundary_t::FLOW,
                        [](const std::array<double, 2>& x) { return 0; },
                    }
                }
            },
            [](const std::array<double, 2>&) { return 0; }, // Правая часть
            bell // Функция влияния
        );

/*
                nonlocal::structural::structural_solver<double, int, int> fem_sol{mesh_proxy};
                nonlocal::structural::equation_parameters<double> parameters;
                parameters.nu = 0.3;
                parameters.E = 21;
                parameters.p1 = p1;
                parameters.type = nonlocal::structural::calc_t::PLANE_STRESS;

                auto sol = fem_sol.stationary(parameters,
                    { // Граничные условия
                        {   "Right",
                            {
                            nonlocal::structural::boundary_t::PRESSURE,
                            [](const std::array<double, 2>& x) { return f1(x); },
                            nonlocal::structural::boundary_t::PRESSURE,
                            [](const std::array<double, 2>&) { return 0; }
                            }
                        },

                        {   "Horizontal",
                                {
                            nonlocal::structural::boundary_t::PRESSURE,
                            [](const std::array<double, 2>&) { return 0; },
                            nonlocal::structural::boundary_t::DISPLACEMENT,
                            [](const std::array<double, 2>&) { return 0; }
                                }
                        },

                        {  "Left",
                                {
                            nonlocal::structural::boundary_t::PRESSURE,
                            [](const std::array<double, 2>& x) { return -f1(x); },
                            nonlocal::structural::boundary_t::PRESSURE,
                            [](const std::array<double, 2>&) { return 0; }
                                }
                        },

                        {   "Vertical",
                            {
                                nonlocal::structural::boundary_t::DISPLACEMENT,
                                [](const std::array<double, 2>&) { return 0; },
                                nonlocal::structural::boundary_t::PRESSURE,
                                [](const std::array<double, 2>&) { return 0; }
                            }
                        }
                    },
                    {}, // Правая часть
                    bell // Функция влияния
                );
*/
                std::cout << std::endl;
                /*
                sol.calc_strain_and_stress();

                if (MPI_utils::MPI_rank() == 0) {
                    std::cout << "Energy = " << sol.calc_energy() << std::endl;
                    using namespace std::string_literals;
                    sol.save_as_vtk(argv[4] + "/structural.vtk"s);
                    save_raw_data(argv[4], *mesh, sol);
                }
                 */
            }
            std::cout << std::endl << std::endl;
            std::cout << std::endl << std::endl;
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