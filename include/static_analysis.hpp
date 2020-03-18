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

enum class boundary_type : uint8_t {TRANSLATION, PRESSURE};

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
                           const std::function<double(double, double)> &right_part,
                           const double p1, const std::function<double(double, double, double, double)> &influence_fun);

std::array<std::vector<double>, 6>
    strains_and_stress(const mesh_2d<double> &mesh, const Eigen::VectorXd &u, const parameters<double> &params,
                       const double p1, const std::function<double(double, double, double, double)> &influence_fun);

void raw_output(const std::string &path,            const mesh_2d<double> &mesh,        const Eigen::VectorXd &u,
                const std::vector<double> &eps11,   const std::vector<double> &eps22,   const std::vector<double> &eps12,
                const std::vector<double> &sigma11, const std::vector<double> &sigma22, const std::vector<double> &sigma12);

void save_as_vtk(const std::string &path,            const mesh_2d<double> &mesh,        const Eigen::VectorXd &u,
                 const std::vector<double> &eps11,   const std::vector<double> &eps22,   const std::vector<double> &eps12,
                 const std::vector<double> &sigma11, const std::vector<double> &sigma22, const std::vector<double> &sigma12);

}