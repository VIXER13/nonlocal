#include <tuple>
#include <functional>
#include "mesh_2d.hpp"
#include "Eigen/Dense"

namespace statics_with_nonloc
{

enum class component : uint8_t {X, Y};

template<class Type>
struct parameters
{
    Type nu, // Коэффициент Пуассона
         E;  // Модуль Юнга
};

enum class boundary_type : uint8_t {DISPLACEMENT, PRESSURE};

template<class Type>
struct boundary_condition
{
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    std::function<Type(Type, Type)> func_x = [](Type, Type) { return 0.; },
                                    func_y = [](Type, Type) { return 0.; };
    boundary_type type_x = boundary_type::PRESSURE,
                  type_y = boundary_type::PRESSURE;
};

Eigen::VectorXd stationary(const mesh_2d<double> &mesh, const parameters<double> &params,
                           const std::vector<boundary_condition<double>> &bounds_cond,
                           const double p1, const std::function<double(double, double, double, double)> &influence_fun);

std::array<std::vector<double>, 6>
    strains_and_stress(const mesh_2d<double> &mesh, const Eigen::VectorXd &u, const parameters<double> &params,
                       const double p1, const std::function<double(double, double, double, double)> &influence_fun);

template<class Type, class Index>
void raw_output(const std::string &path,          const mesh_2d<Type, Index> &mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1> &u,
                const std::vector<Type> &eps11,   const std::vector<Type> &eps22,   const std::vector<Type> &eps12,
                const std::vector<Type> &sigma11, const std::vector<Type> &sigma22, const std::vector<Type> &sigma12)
{
    std::ofstream fout_ux(path + std::string("u_x.csv")),
                  fout_uy(path + std::string("u_y.csv"));
    fout_ux.precision(20);
    fout_uy.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        fout_ux << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i) << std::endl;
        fout_uy << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i+1) << std::endl;
    }

    std::ofstream fout_eps11(path + std::string("eps11.csv")),
                  fout_eps22(path + std::string("eps22.csv")),
                  fout_eps12(path + std::string("eps12.csv"));
    fout_eps11.precision(20);
    fout_eps22.precision(20);
    fout_eps12.precision(20);
    for(size_t i = 0; i < eps11.size(); ++i)
        fout_eps11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps11[i] << std::endl;

    for(size_t i = 0; i < eps22.size(); ++i)
        fout_eps22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps22[i] << std::endl;

    for(size_t i = 0; i < eps12.size(); ++i)
        fout_eps12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps12[i] << std::endl;

    std::ofstream fout_sigma11(path + std::string("sigma11.csv")),
                  fout_sigma22(path + std::string("sigma22.csv")),
                  fout_sigma12(path + std::string("sigma12.csv"));
    fout_sigma11.precision(20);
    fout_sigma22.precision(20);
    fout_sigma12.precision(20);
    for(size_t i = 0; i < sigma11.size(); ++i)
        fout_sigma11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma11[i] << std::endl;

    for(size_t i = 0; i < sigma22.size(); ++i)
        fout_sigma22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma22[i] << std::endl;

    for(size_t i = 0; i < sigma12.size(); ++i)
        fout_sigma12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma12[i] << std::endl;
}

template<class Type, class Index>
void save_as_vtk(const std::string &path,          const mesh_2d<Type, Index> &mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1> &u,
                 const std::vector<Type> &eps11,   const std::vector<Type> &eps22,   const std::vector<Type> &eps12,
                 const std::vector<Type> &sigma11, const std::vector<Type> &sigma22, const std::vector<Type> &sigma12)
{
    std::ofstream fout(path);
    fout.precision(20);

    fout << "# vtk DataFile Version 4.2" << std::endl
         << "Temperature"                << std::endl
         << "ASCII"                      << std::endl
         << "DATASET UNSTRUCTURED_GRID"  << std::endl;

    fout << "POINTS " << mesh.nodes_count() << " double" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << mesh.coord(i, 0) << " " << mesh.coord(i, 1) << " 0" << std::endl;

    fout << "CELLS " << mesh.elements_count() << " " << mesh.elements_count() * 5 << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        fout << 4 << " " << mesh.node_number(i, 0) << " "
                         << mesh.node_number(i, 1) << " "
                         << mesh.node_number(i, 2) << " "
                         << mesh.node_number(i, 3) << std::endl;

    fout << "CELL_TYPES " << mesh.elements_count() << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        fout << 9 << std::endl;

    fout << "POINT_DATA " << mesh.nodes_count() << std::endl;

    fout << "SCALARS U_X double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << u[2*i] << std::endl;

    fout << "SCALARS U_Y double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << u[2*i+1] << std::endl;

    fout << "SCALARS EPS_XX double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps11.size(); ++i)
        fout << eps11[i] << std::endl;

    fout << "SCALARS EPS_YY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps22.size(); ++i)
        fout << eps22[i] << std::endl;

    fout << "SCALARS EPS_XY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps12.size(); ++i)
        fout << eps12[i] << std::endl;

    fout << "SCALARS SIGMA_XX double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma11.size(); ++i)
        fout << sigma11[i] << std::endl;

    fout << "SCALARS SIGMA_YY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma22.size(); ++i)
        fout << sigma22[i] << std::endl;

    fout << "SCALARS SIGMA_XY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma12.size(); ++i)
        fout << sigma12[i] << std::endl;
}

}