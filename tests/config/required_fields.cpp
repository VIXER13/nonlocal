#include "mesh_data.hpp"
#include "time_data.hpp"
#include "boundaries_conditions_data.hpp"

#include "required_fields_json.h"

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::config;

template<class T, size_t N>
class mock_data final {
    nlohmann::json _value;

public:
    explicit constexpr mock_data() noexcept = default;
    explicit constexpr mock_data(const nlohmann::json& value) 
        : _value{value} {}

    constexpr operator nlohmann::json() const {
        return _value;
    }
};

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

template<class T>
auto time_test(const nlohmann::json& config) {
    return [&config] {
        const std::string type_name = '_' + std::string{reflection::type_name<T>()};
        test("time_missed_all" + type_name) = expect_throw<time_data<T>>(config["time_missed_all"]);
        test("time_missed_time_step" + type_name) = expect_throw<time_data<T>>(config["time_missed_time_step"]);
        test("time_missed_steps_count" + type_name) = expect_throw<time_data<T>>(config["time_missed_steps_count"]);
        test("time_all_required_exists" + type_name) = expect_no_throw<time_data<T>>(config["time_all_required_exists"]);
    };
}

template<class T>
auto boundaries_conditions_test(const nlohmann::json& config) {
    return [&config] {
        const std::string type_name = '_' + std::string{reflection::type_name<T>()};
        test("boundaries_conditions_1d_all_missed" + type_name) = expect_throw<boundaries_conditions_data<mock_data, T, 1>>(config["boundaries_conditions_all_missed"]);
        test("boundaries_conditions_1d_missed_left" + type_name) = expect_throw<boundaries_conditions_data<mock_data, T, 1>>(config["boundaries_conditions_missed_left"]);
        test("boundaries_conditions_1d_missed_right" + type_name) = expect_throw<boundaries_conditions_data<mock_data, T, 1>>(config["boundaries_conditions_missed_right"]);
        test("boundaries_conditions_1d_all_required_exists" + type_name) = expect_no_throw<boundaries_conditions_data<mock_data, T, 1>>(config["boundaries_conditions_all_required_exists"]);
        test("boundaries_conditions_2d_all_missed" + type_name) = expect_no_throw<boundaries_conditions_data<mock_data, T, 2>>(config["boundaries_conditions_all_missed"]);
        test("boundaries_conditions_3d_all_missed" + type_name) = expect_no_throw<boundaries_conditions_data<mock_data, T, 3>>(config["boundaries_conditions_all_missed"]);
    };
}

const suite<"config_required_fields"> _ = [] {
    const nlohmann::json config = nlohmann::json::parse(required_fields_json_data);
    
    test("mesh_1d_data_missed_all") = expect_no_throw<mesh_data<1>>(config["mesh_dim_all_required_exists"]);
    test("mesh_2d_data_missed_all") = expect_throw<mesh_data<2>>(config["mesh_dim_missed_all"]);
    test("mesh_2d_data_all_required_exists") = expect_no_throw<mesh_data<2>>(config["mesh_dim_all_required_exists"]);
    test("mesh_3d_data_missed_all") = expect_throw<mesh_data<3>>(config["mesh_dim_missed_all"]);
    test("mesh_3d_data_all_required_exists") = expect_no_throw<mesh_data<3>>(config["mesh_dim_all_required_exists"]);

    test("time_data_float") = time_test<float>(config);
    test("time_data_double") = time_test<double>(config);
    test("time_data_long_double") = time_test<long double>(config);

    test("boundaries_conditions_data_float") = boundaries_conditions_test<float>(config);
    test("boundaries_conditions_data_double") = boundaries_conditions_test<double>(config);
    test("boundaries_conditions_data_long_double") = boundaries_conditions_test<long double>(config);
};

}