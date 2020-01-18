#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <tuple>
#include <functional>
#include "mesh_2d.hpp"

namespace heat_equation_with_nonloc
{

enum class boundary_type {TEMPERATURE, FLOW};

void stationary(const std::string &path, const mesh_2d<double> &mesh,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &boundary,
                const std::function<double(double, double)> &right_part,
                const double p1, const std::function<double(double, double, double, double)> &influence_fun,
                const double volume = 0.0); // В случае задачи Неймана

void nonstationary(const std::string &path,
                   const mesh_2d<double> &mesh, const double tau, const size_t time_steps,
                   const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &boundary,
                   const std::function<double(double, double)> &init_dist,
                   const std::function<double(double, double)> &right_part,
                   const double p1, const std::function<double(double, double, double, double)> &influence_fun,
                   const uint64_t print_frequency = uint64_t(-1));

}

#endif