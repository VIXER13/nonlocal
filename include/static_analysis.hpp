#include <tuple>
#include <functional>
#include "mesh_2d.hpp"
#include "Eigen/Dense"

namespace statics_with_nonloc
{

template<class Type>
struct parameters
{
    Type nu, // Коэффициент Пуассона
         E;  // Модуль Юнга
};

enum class boundary_type {TRANSLATION, FORCE};

Eigen::VectorXd stationary(const std::string &path, const mesh_2d<double> &mesh, const parameters<double> &params,
                           const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                        boundary_type, std::function<double(double, double)>>> &bounds_cond,
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