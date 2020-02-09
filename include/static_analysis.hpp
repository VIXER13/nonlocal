#include <tuple>
#include <functional>
#include "mesh_2d.hpp"

namespace statics_with_nonloc
{

template<class Type>
struct parameters
{
    Type nu, // Коэффициент Пуассона
         E;  // Модуль Юнга
};

enum class boundary_type {TRANSLATION, FORCE};

void stationary(const std::string &path, const mesh_2d<double> &mesh, const parameters<double> &params,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                             boundary_type, std::function<double(double, double)>>> &bounds_cond,
                const double p1, const std::function<double(double, double, double, double)> &influence_fun);
}