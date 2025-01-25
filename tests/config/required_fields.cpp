#include "tests_config_utils.hpp"
#include "required_fields_json.h"

#include "nonlocal_config.hpp"

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace unit_tests;
using namespace nonlocal::config;

template<class T>
auto expect_throw(const nlohmann::json& name) {
    return [&name] {
        expect(throws([&name] { T{name}; }));
    };
}

template<class T>
auto expect_no_throw(const nlohmann::json& name) {
    return [&name] {
        expect(nothrow([&name] { T{name}; }));
    };
}

auto mesh_data_test(const nlohmann::json& config) {
    return [&config] {
        test("mesh_2d_data_missed_all")          = expect_throw<mesh_data<2>>(config["mesh_dim_missed_all"]);
        test("mesh_3d_data_missed_all")          = expect_throw<mesh_data<3>>(config["mesh_dim_missed_all"]);
        test("mesh_1d_data_missed_all")          = expect_no_throw<mesh_data<1>>(config["mesh_dim_all_required_exists"]);
        test("mesh_2d_data_all_required_exists") = expect_no_throw<mesh_data<2>>(config["mesh_dim_all_required_exists"]);
        test("mesh_3d_data_all_required_exists") = expect_no_throw<mesh_data<3>>(config["mesh_dim_all_required_exists"]);
    };
}

template<std::floating_point T>
auto time_data_test(const nlohmann::json& config) {
    return [&config] {
        const std::string suffix = '_' + std::string{reflection::type_name<T>()};
        test("time_missed_all" + suffix)          = expect_throw<time_data<T>>(config["time_missed_all"]);
        test("time_missed_time_step" + suffix)    = expect_throw<time_data<T>>(config["time_missed_time_step"]);
        test("time_missed_steps_count" + suffix)  = expect_throw<time_data<T>>(config["time_missed_steps_count"]);
        test("time_all_required_exists" + suffix) = expect_no_throw<time_data<T>>(config["time_all_required_exists"]);
    };
}

template<std::floating_point T>
auto model_data_test(const nlohmann::json& config) {
    return [&config] {
        const std::string suffix = '_' + std::string{reflection::type_name<T>()};
        test("model_1d_missed_all" + suffix)             = expect_throw<model_data<T, 1>>(config["model_dim_missed_all"]);
        test("model_2d_missed_all" + suffix)             = expect_throw<model_data<T, 2>>(config["model_dim_missed_all"]);
        test("model_3d_missed_all" + suffix)             = expect_throw<model_data<T, 3>>(config["model_dim_missed_all"]);
        test("model_1d_missed_local_weight" + suffix)    = expect_throw<model_data<T, 1>>(config["model_dim_missed_local_weight"]);
        test("model_1d_missed_nonlocal_radius" + suffix) = expect_throw<model_data<T, 1>>(config["model_dim_missed_nonlocal_radius"]);
        test("model_2d_missed_local_weight" + suffix)    = expect_throw<model_data<T, 2>>(config["model_dim_missed_local_weight"]);
        test("model_2d_missed_nonlocal_radius" + suffix) = expect_throw<model_data<T, 2>>(config["model_dim_missed_nonlocal_radius"]);
        test("model_3d_missed_local_weight" + suffix)    = expect_throw<model_data<T, 3>>(config["model_dim_missed_local_weight"]);
        test("model_3d_missed_nonlocal_radius" + suffix) = expect_throw<model_data<T, 3>>(config["model_dim_missed_nonlocal_radius"]);
        test("model_1d_all_required_exists" + suffix)    = expect_no_throw<model_data<T, 1>>(config["model_1d_all_required_exists"]);
        test("model_2d_all_required_exists" + suffix)    = expect_no_throw<model_data<T, 2>>(config["model_2d_all_required_exists"]);
        test("model_3d_all_required_exists" + suffix)    = expect_no_throw<model_data<T, 3>>(config["model_3d_all_required_exists"]);
    };
}

template<std::floating_point T>
auto material_data_test(const nlohmann::json& config) {
    return [&config] {
        const std::string suffix = '_' + std::string{reflection::type_name<T>()};
        test("material_1d_missed_all" + suffix)            = expect_throw<material_data<mock_data, T, 1>>(config["material_dim_missed_all"]);
        test("material_2d_missed_all" + suffix)            = expect_throw<material_data<mock_data, T, 2>>(config["material_dim_missed_all"]);
        test("material_3d_missed_all" + suffix)            = expect_throw<material_data<mock_data, T, 3>>(config["material_dim_missed_all"]);
        test("material_1d_missed_elements_count" + suffix) = expect_throw<material_data<mock_data, T, 1>>(config["material_1d_missed_elements_count"]);
        test("material_1d_missed_length" + suffix)         = expect_throw<material_data<mock_data, T, 1>>(config["material_1d_missed_length"]);
        test("material_1d_missed_physical" + suffix)       = expect_throw<material_data<mock_data, T, 1>>(config["material_1d_missed_physical"]);
        test("material_1d_all_required_exists" + suffix)   = expect_no_throw<material_data<mock_data, T, 1>>(config["material_1d_all_required_exists"]);
        test("material_2d_all_required_exists" + suffix)   = expect_no_throw<material_data<mock_data, T, 2>>(config["material_dim_all_required_exists"]);
        test("material_3d_all_required_exists" + suffix)   = expect_no_throw<material_data<mock_data, T, 3>>(config["material_dim_all_required_exists"]);
    };
}

template<std::floating_point T>
auto mechanical_material_data_test(const nlohmann::json& config) {
    return [&config] {
        const std::string suffix = '_' + std::string{reflection::type_name<T>()};
        test("mechanical_material_1d_missed_all" + suffix)                                 = expect_throw<mechanical_material_data<T, 1>>(config["mechanical_material_dim_missed_all"]);
        test("mechanical_material_1d_all_required_exists" + suffix)                        = expect_no_throw<mechanical_material_data<T, 1>>(config["mechanical_material_1d_all_required_exists"]);
        test("mechanical_material_2d_isotropic_all_params" + suffix)                       = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_isotropic_all_params"]);
        test("mechanical_material_2d_isotropic_all_required_exists" + suffix)              = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_isotropic_all_required_exists"]);
        test("mechanical_material_2d_orthotropic_all_params_exists" + suffix)              = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_all_params_exists"]);
        test("mechanical_material_2d_orthotropic_two_nulls_first" + suffix)                = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_first"]);
        test("mechanical_material_2d_orthotropic_two_nulls_second" + suffix)               = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_second"]);
        test("mechanical_material_2d_orthotropic_two_nulls_third" + suffix)                = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_third"]);
        test("mechanical_material_2d_orthotropic_two_nulls_fourth" + suffix)               = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_fourth"]);
        test("mechanical_material_2d_orthotropic_two_nulls_fifth" + suffix)                = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_fifth"]);
        test("mechanical_material_2d_orthotropic_two_nulls_sixth" + suffix)                = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_two_nulls_sixth"]);
        test("mechanical_material_2d_orthotropic_without_shear_modulus_first" + suffix)    = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_without_shear_modulus_first"]);
        test("mechanical_material_2d_orthotropic_without_shear_modulus_second" + suffix)   = expect_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_without_shear_modulus_second"]);
        test("mechanical_material_2d_orthotropic_youngs_modulus" + suffix)                 = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_youngs_modulus"]);
        test("mechanical_material_2d_orthotropic_poissons_ratio" + suffix)                 = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_poissons_ratio"]);
        test("mechanical_material_2d_orthotropic_one_null_youngs_modulus_first" + suffix)  = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_one_null_youngs_modulus_first"]);
        test("mechanical_material_2d_orthotropic_one_null_youngs_modulus_second" + suffix) = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_one_null_youngs_modulus_second"]);
        test("mechanical_material_2d_orthotropic_one_null_poissons_ratio_first" + suffix)  = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_one_null_poissons_ratio_first"]);
        test("mechanical_material_2d_orthotropic_one_null_poissons_ratio_second" + suffix) = expect_no_throw<mechanical_material_data<T, 2>>(config["mechanical_material_2d_orthotropic_one_null_poissons_ratio_second"]);
    };
}

template<std::floating_point T>
auto thermal_material_data_test(const nlohmann::json& config) {
    return [&config] {
        const std::string suffix = '_' + std::string{reflection::type_name<T>()};
        test("thermal_material_1d_missed_all" + suffix)                      = expect_throw<thermal_material_data<T, 1>>(config["thermal_material_dim_missed_all"]);
        test("thermal_material_2d_missed_all" + suffix)                      = expect_throw<thermal_material_data<T, 2>>(config["thermal_material_dim_missed_all"]);
        test("thermal_material_1d_all_required_exists" + suffix)             = expect_no_throw<thermal_material_data<T, 1>>(config["thermal_material_dim_all_required_exists"]);
        test("thermal_material_2d_all_required_exists" + suffix)             = expect_no_throw<thermal_material_data<T, 2>>(config["thermal_material_dim_all_required_exists"]);
        test("thermal_material_2d_orthotropic_all_required_exists" + suffix) = expect_no_throw<thermal_material_data<T, 2>>(config["thermal_material_2d_orthotropic_all_required_exists"]);
        test("thermal_material_2d_anisotropic_all_required_exists" + suffix) = expect_no_throw<thermal_material_data<T, 2>>(config["thermal_material_2d_anisotropic_all_required_exists"]);
    };
}

const suite<"config_required_fields"> _ = [] {
    const nlohmann::json config = nlohmann::json::parse(required_fields_json_data);
    test("mesh_data") = mesh_data_test(config);
    test("time_data") = time_data_test<double>(config);
    test("model_data") = model_data_test<double>(config);
    test("material_data") = material_data_test<double>(config);
    test("mechanical_material") = mechanical_material_data_test<double>(config);
    test("thermal_material") = thermal_material_data_test<double>(config);
};

}