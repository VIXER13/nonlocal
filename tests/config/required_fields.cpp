#include "required_fields_json.h"

#include "thermal_config_data.hpp"

#include <gtest/gtest.h>

namespace {

using namespace nonlocal::config;

using model_data_1d = model_data<float, 1u>;
using model_data_2d = model_data<float, 2u>;

template<class T, size_t N>
struct physic_data_mock final {
    explicit constexpr physic_data_mock() noexcept = default;
    explicit physic_data_mock(const Json::Value&) noexcept {}
};

using material_data_1d = material_data<physic_data_mock, float, 1u>;
using material_data_2d = material_data<physic_data_mock, float, 2u>;

using mesh_data_2d = mesh_data<2u>;

const Json::Value tests = read_json(required_fields_json_data, std::next(required_fields_json_data, required_fields_json_size));

}



TEST(config_required_fields, nonstationary_missed_all) {
    EXPECT_THROW(nonstationary_data<float>{tests["nonstationary_missed_all"]}, std::exception);
}

TEST(config_required_fields, nonstationary_missed_time_step) {
    EXPECT_THROW(nonstationary_data<float>{tests["nonstationary_missed_time_step"]}, std::exception);
}

TEST(config_required_fields, nonstationary_missed_steps_count) {
    EXPECT_THROW(nonstationary_data<float>{tests["nonstationary_missed_steps_count"]}, std::exception);
}

TEST(config_required_fields, nonstationary_all_required_exists) {
    EXPECT_NO_THROW(nonstationary_data<float>{tests["nonstationary_all_required_exists"]});
}



TEST(config_required_fields, model_data_1d_missed_all) {
    EXPECT_THROW(model_data_1d{tests["model_data_1d_missed_all"]}, std::exception);
}

TEST(config_required_fields, model_data_1d_missed_local_weight) {
    EXPECT_THROW(model_data_1d{tests["model_data_1d_missed_local_weight"]}, std::exception);
}

TEST(config_required_fields, model_data_1d_missed_nonlocal_radius) {
    EXPECT_THROW(model_data_1d{tests["model_data_1d_missed_nonlocal_radius"]}, std::exception);
}

TEST(config_required_fields, model_data_1d_all_required_exists)
{
    EXPECT_NO_THROW(model_data_1d{tests["model_data_1d_all_required_exists"]});
}



TEST(config_required_fields, model_data_2d_missed_all) {
    EXPECT_THROW(model_data_2d{tests["model_data_2d_missed_all"]}, std::exception);
}

TEST(config_required_fields, model_data_2d_missed_local_weight) {
    EXPECT_THROW(model_data_2d{tests["model_data_2d_missed_local_weight"]}, std::exception);
}

TEST(config_required_fields, model_data_2d_missed_nonlocal_radius) {
    EXPECT_THROW(model_data_2d{tests["model_data_2d_missed_nonlocal_radius"]}, std::exception);
}

TEST(config_required_fields, model_data_2d_all_required_exists) {
    EXPECT_NO_THROW(model_data_2d{tests["model_data_2d_all_required_exists"]});
}



TEST(config_required_fields, material_data_1d_missed_all) {
    EXPECT_THROW(material_data_1d{tests["material_data_1d_missed_all"]}, std::exception);
}

TEST(config_required_fields, material_data_1d_missed_elements_count) {
    EXPECT_THROW(material_data_1d{tests["material_data_1d_missed_elements_count"]}, std::exception);
}

TEST(config_required_fields, material_data_1d_missed_length) {
    EXPECT_THROW(material_data_1d{tests["material_data_1d_missed_length"]}, std::exception);
}

TEST(config_required_fields, material_data_1d_missed_physical) {
    EXPECT_THROW(material_data_1d{tests["material_data_1d_missed_physical"]}, std::exception);
}

TEST(config_required_fields, material_data_1d_all_required_exists) {
    EXPECT_NO_THROW(material_data_1d{tests["material_data_1d_all_required_exists"]});
}



TEST(config_required_fields, material_data_2d_missed_all) {
    EXPECT_THROW(material_data_2d{tests["material_data_2d_missed_all"]}, std::exception);
}

TEST(config_required_fields, material_data_2d_all_required_exists) {
    EXPECT_NO_THROW(material_data_2d{tests["material_data_2d_all_required_exists"]});
}



TEST(config_required_fields, mesh_data_2d_missed_all) {
    EXPECT_THROW(mesh_data_2d{tests["mesh_data_2d_missed_all"]}, std::exception);
}

TEST(config_required_fields, mesh_data_2d_all_required_exists) {
    EXPECT_NO_THROW(mesh_data_2d{tests["mesh_data_2d_all_required_exists"]});
}