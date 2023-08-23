#ifndef NONLOCFEM_INIT_UTILS_HPP
#define NONLOCFEM_INIT_UTILS_HPP

#include "nonlocal_config.hpp"

#include <set>
#include <iostream>
#include <exception>

namespace nonlocal {

void init_save_data(const nonlocal::config::save_data& save, const nlohmann::json& config);
std::string init_available_problems_list(const std::set<config::problem_t>& available_problems);

bool is_thermal_problem(const config::problem_t problem);
bool is_mechanical_problem(const config::problem_t problem);

}

#endif