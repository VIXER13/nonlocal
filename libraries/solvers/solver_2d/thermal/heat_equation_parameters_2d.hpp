#ifndef HEAT_EQUATION_PARAMETERS_2D_HPP
#define HEAT_EQUATION_PARAMETERS_2D_HPP

#include "../../solvers_constants.hpp"
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace nonlocal::thermal {

template<class T, material_t Material>
struct equation_parameters_2d final {
    static_assert(Material == material_t::ISOTROPIC ||
                  Material == material_t::ORTHOTROPIC,
                  "Only isotropic and orthotropic materials are supported");

    std::conditional_t<
        Material == material_t::ISOTROPIC,
        T,
        std::array<T, 2>
    > thermal_conductivity;
    std::unordered_map<std::string, T> heat_transfer;
    T heat_capacity = T{1};
    T density = T{1};
    T integral = T{0}; // Solution integral for Neumann problem

    explicit equation_parameters_2d(const std::vector<std::string>& names = {}) noexcept {
        if constexpr(Material == material_t::ISOTROPIC)
            thermal_conductivity = T{1};
        if constexpr(Material == material_t::ORTHOTROPIC)
            thermal_conductivity.fill(T{1});
        for(const std::string& name : names)
            heat_transfer[name] = T{1};
    }
};

}

#endif