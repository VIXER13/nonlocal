#include "tests_config_utils.hpp"
#include "double_conversion_json.h"

#include "nonlocal_config.hpp"

#include <boost/ut.hpp>

#include <iostream>

namespace {

using namespace boost::ut;
using namespace unit_tests;
using namespace nonlocal::config;

template<class T>
auto double_conversion(const nlohmann::json& config) {
    return [&config] {
        const T data{config};
        const nlohmann::json back_converted = data;
        expect(eq(config, back_converted));
    };
}
    
const suite<"config_double_conversion"> _ = [] {
    const nlohmann::json config = nlohmann::json::parse(double_conversion_json_data);
    test("mesh_2d") = double_conversion<mesh_data<2>>(config["mesh_2d"]);
    test("time") = double_conversion<time_data<double>>(config["time"]);
    test("model_1d") = double_conversion<model_data<double, 1>>(config["model_1d"]);
    test("model_2d") = double_conversion<model_data<double, 2>>(config["model_2d"]);
    test("material_1d") = double_conversion<material_data<mock_data, double, 1>>(config["material_1d"]);
    test("material_2d") = double_conversion<material_data<mock_data, double, 2>>(config["material_2d"]);
    test("thermal_auxiliary") = double_conversion<thermal_auxiliary_data<double>>(config["thermal_auxiliary"]);
    test("thermal_material_1d") = double_conversion<thermal_material_data<double, 1>>(config["thermal_material_1d"]);
    for(const std::string material : {"thermal_isotropic_material_2d", "thermal_orthotropic_material_2d", "thermal_anisotropic_material_2d"})
        test(material) = double_conversion<thermal_material_data<double, 2>>(config[material]);
};

}