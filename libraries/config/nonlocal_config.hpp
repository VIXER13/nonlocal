#ifndef NONLOCAL_CONFIG_HPP
#define NONLOCAL_CONFIG_HPP

#include "task_data.hpp"
#include "save_data.hpp"
#include "mesh_data.hpp"
#include "time_data.hpp"
#include "thermal_auxiliary_data.hpp"
#include "boundaries_conditions_data.hpp"
#include "thermal_boundary_condition_data.hpp"
#include "mechanical_boundary_condition_data.hpp"
#include "materials_data.hpp"
#include "thermal_material_data.hpp"
#include "mechanical_material_data.hpp"

namespace nonlocal::config {

template<class T>
using thermal_boundaries_conditions_1d = boundaries_conditions_data<thermal_boundary_condition_data, T, 1>;

template<class T>
using thermal_boundaries_conditions_2d = boundaries_conditions_data<thermal_boundary_condition_data, T, 2>;

template<class T>
using mechanical_boundaries_conditions_1d = boundaries_conditions_data<mechanical_boundary_condition_data, T, 1>;

template<class T>
using mechanical_boundaries_conditions_2d = boundaries_conditions_data<mechanical_boundary_condition_data, T, 2>;

template<class T>
using thermal_materials_1d = materials_data<thermal_material_data, T, 1>;

template<class T>
using thermal_materials_2d = materials_data<thermal_material_data, T, 2>;

template<class T>
using mechanical_materials_1d = materials_data<mechanical_material_data, T, 1>;

template<class T>
using mechanical_materials_2d = materials_data<mechanical_material_data, T, 2>;

}

#endif