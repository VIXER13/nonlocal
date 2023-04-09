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

template<class Data>
bool exception_catched(const Json::Value& value) {
    try {
        const Data data{value};
    } catch(const std::exception&) {
        return true;
    }
    return false;
}

}



TEST(config_required_fields, nonstationary_missed_all)
{
    EXPECT_TRUE(exception_catched<nonstationary_data<float>>(tests["nonstationary_missed_all"]));
}

TEST(config_required_fields, nonstationary_missed_time_step)
{
    EXPECT_TRUE(exception_catched<nonstationary_data<float>>(tests["nonstationary_missed_time_step"]));
}

TEST(config_required_fields, nonstationary_missed_steps_count)
{
    EXPECT_TRUE(exception_catched<nonstationary_data<float>>(tests["nonstationary_missed_steps_count"]));
}

TEST(config_required_fields, nonstationary_all_required_exists)
{
    EXPECT_FALSE(exception_catched<nonstationary_data<float>>(tests["nonstationary_all_required_exists"]));
}



TEST(config_required_fields, model_data_1d_missed_all)
{
    EXPECT_TRUE(exception_catched<model_data_1d>(tests["model_data_1d_missed_all"]));
}

TEST(config_required_fields, model_data_1d_missed_local_weight)
{
    EXPECT_TRUE(exception_catched<model_data_1d>(tests["model_data_1d_missed_local_weight"]));
}

TEST(config_required_fields, model_data_1d_missed_nonlocal_radius)
{
    EXPECT_TRUE(exception_catched<model_data_1d>(tests["model_data_1d_missed_nonlocal_radius"]));
}

TEST(config_required_fields, model_data_1d_all_required_exists)
{
    EXPECT_FALSE(exception_catched<model_data_1d>(tests["model_data_1d_all_required_exists"]));
}



TEST(config_required_fields, model_data_2d_missed_all)
{
    EXPECT_TRUE(exception_catched<model_data_2d>(tests["model_data_2d_missed_all"]));
}

TEST(config_required_fields, model_data_2d_missed_local_weight)
{
    EXPECT_TRUE(exception_catched<model_data_2d>(tests["model_data_2d_missed_local_weight"]));
}

TEST(config_required_fields, model_data_2d_missed_nonlocal_radius)
{
    EXPECT_TRUE(exception_catched<model_data_2d>(tests["model_data_2d_missed_nonlocal_radius"]));
}

TEST(config_required_fields, model_data_2d_all_required_exists)
{
    EXPECT_FALSE(exception_catched<model_data_2d>(tests["model_data_2d_all_required_exists"]));
}



TEST(config_required_fields, material_data_1d_missed_all)
{
    EXPECT_TRUE(exception_catched<material_data_1d>(tests["material_data_1d_missed_all"]));
}

TEST(config_required_fields, material_data_1d_missed_elements_count)
{
    EXPECT_TRUE(exception_catched<material_data_1d>(tests["material_data_1d_missed_elements_count"]));
}

TEST(config_required_fields, material_data_1d_missed_length)
{
    EXPECT_TRUE(exception_catched<material_data_1d>(tests["material_data_1d_missed_length"]));
}

TEST(config_required_fields, material_data_1d_missed_physical)
{
    EXPECT_TRUE(exception_catched<material_data_1d>(tests["material_data_1d_missed_physical"]));
}

TEST(config_required_fields, material_data_1d_all_required_exists)
{
    EXPECT_FALSE(exception_catched<material_data_1d>(tests["material_data_1d_all_required_exists"]));
}



TEST(config_required_fields, material_data_2d_missed_all)
{
    EXPECT_TRUE(exception_catched<material_data_2d>(tests["material_data_2d_missed_all"]));
}

TEST(config_required_fields, material_data_2d_all_required_exists)
{
    EXPECT_FALSE(exception_catched<material_data_2d>(tests["material_data_2d_all_required_exists"]));
}



TEST(config_required_fields, mesh_data_2d_missed_all)
{
    EXPECT_TRUE(exception_catched<mesh_data_2d>(tests["mesh_data_2d_missed_all"]));
}

TEST(config_required_fields, mesh_data_2d_all_required_exists)
{
    EXPECT_FALSE(exception_catched<mesh_data_2d>(tests["mesh_data_2d_all_required_exists"]));
}