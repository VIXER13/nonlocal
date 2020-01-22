#include <tuple>
#include <functional>
#include "mesh_2d.hpp"

enum class boundary_type {FIXED, FORCE};

void stationary(const std::string &path, const mesh_2d<double> &mesh,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                             boundary_type, std::function<double(double, double)>>> &bounds_cond);