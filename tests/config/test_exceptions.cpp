#include "test_exceptions_json.h"

#include <config/read_mechanical_boundary_conditions.hpp>
#include <config/read_thermal_boundary_conditions.hpp>
#include <config/read_mesh.hpp>
#include <config/read_model.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::config;

template<class Function>
void expect_throws(const Function& function, const nlohmann::json& config, const std::string& field) {
    expect(throws<std::domain_error>([&function, &config, &field]{ function(config[field], field); }));
}

template<class Function>
void expect_nothrows(const Function& function, const nlohmann::json& config, const std::string& field) {
    expect(nothrow([&function, &config, &field]{ function(config[field], field); }));
}

const suite<"config_exceptions"> _ = [] {
    using T = double;
    const nlohmann::json config = nlohmann::json::parse(test_exceptions_json_data);

    "read_model_1d"_test = [&config]{
        expect_throws(read_model_1d<T>, config, "model_1d_local_weight_is_missed_fail");
        expect_throws(read_model_1d<T>, config, "model_1d_nonlocal_radius_is_missed_fail");
        expect_throws(read_model_1d<T>, config, "model_1d_local_weight_is_negative_fail");
        expect_throws(read_model_1d<T>, config, "model_1d_local_weight_is_zero_fail");
        expect_throws(read_model_1d<T>, config, "model_1d_local_weight_greater_1_fail");
        expect_throws(read_model_1d<T>, config, "model_1d_nonlocal_radius_is_negative_fail");
        expect_nothrows(read_model_1d<T>, config, "model_1d_ok");
    };

    "read_model_2d"_test = [&config]{
        expect_throws(read_model_2d<T>, config, "model_2d_local_weight_is_missed_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_nonlocal_radius_is_missed_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_nonlocal_radius_wrong_dimension_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_local_weight_is_negative_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_local_weight_is_zero_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_local_weight_greater_1_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_nonlocal_radius_is_negative_fail");
        expect_throws(read_model_2d<T>, config, "model_2d_nonlocal_one_of_the_radius_is_negative_fail");
        expect_nothrows(read_model_2d<T>, config, "model_2d_ok");
        expect_nothrows(read_model_2d<T>, config, "model_2d_two_radii_ok");
    };

    "read_thermal_boundaries_conditions_1d"_test = [&config]{
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_empty_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_left_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_right_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_kind_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_flux_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_third_kind_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_third_kind_missed_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_third_kind_negative_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_missed_emissivity_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_negative_emissivity_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_emissivity_greater_than_1_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_combined_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_combined_negative_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_combined_negative_emissivity_fail");
        expect_nothrows(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_temperature_and_flux_ok");
        expect_nothrows(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_convenction_and_radiation_ok");
        expect_nothrows(read_thermal_boundaries_conditions_1d<T>, config, "thermal_boundaries_conditions_combined_ok");
    };

    "read_thermal_boundaries_conditions_2d"_test = [&config]{
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_empty_fail");
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_left_fail");
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_right_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_kind_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_flux_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_third_kind_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_third_kind_missed_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_third_kind_negative_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_missed_emissivity_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_negative_emissivity_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_emissivity_greater_than_1_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_combined_missed_temperature_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_combined_negative_heat_transfer_fail");
        expect_throws(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_combined_negative_emissivity_fail");
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_temperature_and_flux_ok");
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_convenction_and_radiation_ok");
        expect_nothrows(read_thermal_boundaries_conditions_2d<T>, config, "thermal_boundaries_conditions_combined_ok");
    };

    "read_mechanical_boundaries_conditions_2d"_test = [&config]{
        expect_nothrows(read_mechanical_boundaries_conditions_2d<T>, config, "mechanical_boundaries_conditions_2d_empty_ok");
        expect_throws(read_mechanical_boundaries_conditions_2d<T>, config, "mechanical_boundaries_conditions_2d_empty_condition_fail");
        expect_throws(read_mechanical_boundaries_conditions_2d<T>, config, "mechanical_boundaries_conditions_2d_wrong_dimension_1_fail");
        expect_throws(read_mechanical_boundaries_conditions_2d<T>, config, "mechanical_boundaries_conditions_2d_wrong_dimension_3_fail");
        expect_nothrows(read_mechanical_boundaries_conditions_2d<T>, config, "mechanical_boundaries_conditions_2d_ok");
    };

    "read_mesh_1d"_test = [&config] {
        expect_throws(read_mesh_1d<T>, config, "read_mesh_1d_missed_materials_fail");
        expect_throws(read_mesh_1d<T>, config, "read_mesh_1d_materials_not_array_fail");
        expect_throws(read_mesh_1d<T>, config, "read_mesh_1d_missed_elements_count_fail");
        expect_throws(read_mesh_1d<T>, config, "read_mesh_1d_missed_length_fail");
        expect_nothrows(read_mesh_1d<T>, config, "read_mesh_1d_ok");
    };
};

}