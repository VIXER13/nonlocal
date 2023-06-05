#include "tests_config_utils.hpp"

#include "save_data.hpp"
#include "mesh_data.hpp"
#include "time_data.hpp"
#include "boundaries_conditions_data.hpp"
#include "model_data.hpp"
#include "material_data.hpp"
#include "thermal_boundary_condition_data.hpp"
#include "thermal_equation_data.hpp"
#include "thermal_material_data.hpp"

#include "reverse_conversion_json.h"

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace unit_tests;
using namespace nonlocal::config;

template<class T>
auto reverse_conversion(const nlohmann::json& config) {
    return [&config] {
        const T data{config};
        const nlohmann::json back_converted = data;
        expect(eq(config, back_converted));
    };
}
    
const suite<"config_reverse_conversion"> _ = [] {
    const nlohmann::json config = nlohmann::json::parse(reverse_conversion_json_data);

    test("save") = reverse_conversion<save_data>(config["save"]);
    test("save_with_precision") = reverse_conversion<save_data>(config["save_with_precision"]);

    test("mesh_1d") = reverse_conversion<mesh_data<1>>(config["mesh_1d"]);
    test("mesh_2d") = reverse_conversion<mesh_data<2>>(config["mesh_2d"]);

    test("time_float") = reverse_conversion<time_data<float>>(config["time"]);
    test("time_double") = reverse_conversion<time_data<double>>(config["time"]);
    test("time_long_double") = reverse_conversion<time_data<long double>>(config["time"]);

    test("boundaries_conditions_1d_float") = reverse_conversion<boundaries_conditions_data<mock_data, float, 1>>(config["boundaries_conditions_1d"]);
    test("boundaries_conditions_1d_double") = reverse_conversion<boundaries_conditions_data<mock_data, double, 1>>(config["boundaries_conditions_1d"]);
    test("boundaries_conditions_1d_long_double") = reverse_conversion<boundaries_conditions_data<mock_data, long double, 1>>(config["boundaries_conditions_1d"]);

    test("boundaries_conditions_2d_float") = reverse_conversion<boundaries_conditions_data<mock_data, float, 2>>(config["boundaries_conditions_2d"]);
    test("boundaries_conditions_2d_double") = reverse_conversion<boundaries_conditions_data<mock_data, double, 2>>(config["boundaries_conditions_2d"]);
    test("boundaries_conditions_2d_long_double") = reverse_conversion<boundaries_conditions_data<mock_data, long double, 2>>(config["boundaries_conditions_2d"]);

    test("model_1d_float") = reverse_conversion<model_data<float, 1>>(config["model_1d"]);
    test("model_1d_double") = reverse_conversion<model_data<double, 1>>(config["model_1d"]);
    test("model_1d_long_double") = reverse_conversion<model_data<long double, 1>>(config["model_1d"]);

    test("model_2d_float") = reverse_conversion<model_data<float, 2>>(config["model_2d"]);
    test("model_2d_double") = reverse_conversion<model_data<double, 2>>(config["model_2d"]);
    test("model_2d_long_double") = reverse_conversion<model_data<long double, 2>>(config["model_2d"]);

    test("material_1d_float") = reverse_conversion<material_data<mock_data, float, 1>>(config["material_1d"]);
    test("material_1d_double") = reverse_conversion<material_data<mock_data, double, 1>>(config["material_1d"]);
    test("material_1d_long_double") = reverse_conversion<material_data<mock_data, long double, 1>>(config["material_1d"]);

    test("material_2d_float") = reverse_conversion<material_data<mock_data, float, 2>>(config["material_2d"]);
    test("material_2d_double") = reverse_conversion<material_data<mock_data, double, 2>>(config["material_2d"]);
    test("material_2d_long_double") = reverse_conversion<material_data<mock_data, long double, 2>>(config["material_2d"]);

    test("thermal_boundary_condition_float") = reverse_conversion<thermal_boundary_condition_data<float>>(config["thermal_boundary_condition"]);
    test("thermal_boundary_condition_double") = reverse_conversion<thermal_boundary_condition_data<double>>(config["thermal_boundary_condition"]);
    test("thermal_boundary_condition_long_double") = reverse_conversion<thermal_boundary_condition_data<long double>>(config["thermal_boundary_condition"]);

    test("thermal_equation_float") = reverse_conversion<thermal_equation_data<float>>(config["thermal_equation"]);
    test("thermal_equation_double") = reverse_conversion<thermal_equation_data<double>>(config["thermal_equation"]);
    test("thermal_equation_long_double") = reverse_conversion<thermal_equation_data<long double>>(config["thermal_equation"]);

    test("thermal_material_1d_float") = reverse_conversion<thermal_material_data<float, 1>>(config["thermal_material_1d"]);
    test("thermal_material_1d_double") = reverse_conversion<thermal_material_data<double, 1>>(config["thermal_material_1d"]);
    test("thermal_material_1d_long_double") = reverse_conversion<thermal_material_data<long double, 1>>(config["thermal_material_1d"]);

    for(const std::string material : {"thermal_isotropic_material_2d", "thermal_orthotropic_material_2d", "thermal_anisotropic_material_2d"}) {
        test(material + "_float") = reverse_conversion<thermal_material_data<float, 2>>(config[material]);
        test(material + "_double") = reverse_conversion<thermal_material_data<double, 2>>(config[material]);
        test(material + "_long_double") = reverse_conversion<thermal_material_data<long double, 2>>(config[material]);
    }
};

}