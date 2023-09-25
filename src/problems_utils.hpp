#ifndef NONLOCFEM_PROBLEMS_UTILS_2D_HPP
#define NONLOCFEM_PROBLEMS_UTILS_2D_HPP

#include "nonlocal_config.hpp"

namespace nonlocal {

template<std::floating_point T, template<class, size_t> class Physics>
std::unordered_map<std::string, T> get_search_radii(const config::materials_data<Physics, T, 2>& materials) {
    std::unordered_map<std::string, T> result;
    for(const auto& [name, material] : materials.materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            result[name] = std::max(material.model.search_radius[X], material.model.search_radius[Y]);
    return result;
}

}

#endif