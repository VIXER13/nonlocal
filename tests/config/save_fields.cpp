#include "config_test_types.hpp"

#include "save_fields_json.h"

#include <gtest/gtest.h>

namespace {

using namespace nonlocal::config;
using namespace nonlocal::config::test;

const Json::Value tests = read_json(save_fields_json_data, std::next(save_fields_json_data, save_fields_json_size));

}

#define SAVE_TEST(TYPE, FIELD) \
TEST(config_save_fields, FIELD) { \
    const Json::Value& value = tests[#FIELD]; \
    const TYPE data{value}; \
    EXPECT_EQ(Json::Value{data}, value); \
}

SAVE_TEST(save_data, save);
SAVE_TEST(save_data, save_with_precision);
SAVE_TEST(mesh_data_1d, mesh_1d);
SAVE_TEST(mesh_data_2d, mesh_2d);
SAVE_TEST(time_data<float>, time);
SAVE_TEST(boundaries_conditions_data_1d, boundaries_conditions_1d);
SAVE_TEST(boundaries_conditions_data_2d, boundaries_conditions_2d);
SAVE_TEST(model_data_1d, model_1d);
SAVE_TEST(model_data_2d, model_2d);
SAVE_TEST(material_data_1d, material_1d);
SAVE_TEST(material_data_2d, material_2d);
SAVE_TEST(thermal_boundary_condition_1d, thermal_boundary_condition);
SAVE_TEST(thermal_equation_data<float>, thermal_equation);
SAVE_TEST(thermal_material_data_1d, thermal_material_1d);
SAVE_TEST(thermal_material_data_2d, thermal_isotropic_material_2d);
SAVE_TEST(thermal_material_data_2d, thermal_orthotropic_material_2d);
SAVE_TEST(thermal_material_data_2d, thermal_anisotropic_material_2d);

#undef SAVE_TEST