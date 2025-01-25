#pragma once

#include "mesh_data.hpp"
#include "time_data.hpp"
#include "thermal_auxiliary_data.hpp"
#include "materials_data.hpp"
#include "thermal_material_data.hpp"
#include "mechanical_material_data.hpp"

namespace nonlocal::config {

template<class T>
using thermal_materials_1d = materials_data<thermal_material_data, T, 1>;

template<class T>
using thermal_materials_2d = materials_data<thermal_material_data, T, 2>;

template<class T>
using mechanical_materials_1d = materials_data<mechanical_material_data, T, 1>;

template<class T>
using mechanical_materials_2d = materials_data<mechanical_material_data, T, 2>;

}