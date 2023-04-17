#include "config_test_types.hpp"

#include "required_fields_json.h"

#include <gtest/gtest.h>

namespace {

using namespace nonlocal::config;
using namespace nonlocal::config::test;

const Json::Value tests = read_json(required_fields_json_data, std::next(required_fields_json_data, required_fields_json_size));

}

#define THROW_TEST(TYPE, FIELD, GTEST_MACRO) \
TEST(config_required_fields, FIELD) { \
    GTEST_MACRO(TYPE{tests[#FIELD]}); \
}
#define EXPECT_THROW_TEST(TYPE, FIELD) THROW_TEST(TYPE, FIELD, EXPECT_ANY_THROW)
#define EXPECT_NO_THROW_TEST(TYPE, FIELD) THROW_TEST(TYPE, FIELD, EXPECT_NO_THROW)


EXPECT_THROW_TEST(mesh_data_2d, mesh_2d_missed_all);
EXPECT_NO_THROW_TEST(mesh_data_2d, mesh_2d_all_required_exists);

EXPECT_THROW_TEST(boundaries_conditions_data_1d, boundaries_conditions_1d_all_missed);
EXPECT_THROW_TEST(boundaries_conditions_data_1d, boundaries_conditions_1d_missed_left);
EXPECT_THROW_TEST(boundaries_conditions_data_1d, boundaries_conditions_1d_missed_right);
EXPECT_NO_THROW_TEST(boundaries_conditions_data_1d, boundaries_conditions_1d_all_required_exists);

EXPECT_THROW_TEST(time_data<float>, time_missed_all);
EXPECT_THROW_TEST(time_data<float>, time_missed_time_step);
EXPECT_THROW_TEST(time_data<float>, time_missed_steps_count);
EXPECT_NO_THROW_TEST(time_data<float>, time_all_required_exists);

EXPECT_THROW_TEST(model_data_1d, model_1d_missed_all);
EXPECT_THROW_TEST(model_data_1d, model_1d_missed_local_weight);
EXPECT_THROW_TEST(model_data_1d, model_1d_missed_nonlocal_radius);
EXPECT_NO_THROW_TEST(model_data_1d, model_1d_all_required_exists);

EXPECT_THROW_TEST(model_data_2d, model_2d_missed_all);
EXPECT_THROW_TEST(model_data_2d, model_2d_missed_local_weight);
EXPECT_THROW_TEST(model_data_2d, model_2d_missed_nonlocal_radius);
EXPECT_NO_THROW_TEST(model_data_2d, model_2d_all_required_exists);

EXPECT_THROW_TEST(material_data_1d, material_1d_missed_all);
EXPECT_THROW_TEST(material_data_1d, material_1d_missed_elements_count);
EXPECT_THROW_TEST(material_data_1d, material_1d_missed_length);
EXPECT_THROW_TEST(material_data_1d, material_1d_missed_physical);
EXPECT_NO_THROW_TEST(material_data_1d, material_1d_all_required_exists);

EXPECT_THROW_TEST(material_data_2d, material_2d_missed_all);
EXPECT_NO_THROW_TEST(material_data_2d, material_2d_all_required_exists);

EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, thermal_boundary_condition_missed_all);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, thermal_boundary_condition_missed_kind);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, temperature_condition_missed_temperature);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, flux_condition_missed_flux);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, convection_condition_missed_temperature);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, convection_condition_missed_heat_transfer);
EXPECT_THROW_TEST(thermal_boundary_condition_data<float>, radiation_condition_missed_emissivity);
EXPECT_NO_THROW_TEST(thermal_boundary_condition_data<float>, temperature_condition_all_required_exists);
EXPECT_NO_THROW_TEST(thermal_boundary_condition_data<float>, flux_condition_all_required_exists);
EXPECT_NO_THROW_TEST(thermal_boundary_condition_data<float>, convection_condition_all_required_exists);
EXPECT_NO_THROW_TEST(thermal_boundary_condition_data<float>, radiation_condition_all_required_exists);
EXPECT_NO_THROW_TEST(thermal_boundary_condition_data<float>, combined_condition_flux_missed_all);

#undef THROW_TEST
#undef EXPECT_THROW_TEST
#undef EXPECT_NO_THROW_TEST