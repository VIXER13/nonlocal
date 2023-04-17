#include "thermal_config_data.hpp"

namespace nonlocal::config::test {

template<class T, size_t N>
class data_mock final {
    Json::Value _value;

public:
    explicit constexpr data_mock() noexcept = default;
    explicit constexpr data_mock(const Json::Value& value) 
        : _value{value} {}

    constexpr operator Json::Value() const {
        return _value;
    }
};

using mesh_data_1d = mesh_data<1u>;
using mesh_data_2d = mesh_data<2u>;

using model_data_1d = model_data<float, 1u>;
using model_data_2d = model_data<float, 2u>;

using boundaries_conditions_data_1d = boundaries_conditions_data<data_mock, float, 1u>;
using boundaries_conditions_data_2d = boundaries_conditions_data<data_mock, float, 2u>;

using thermal_boundary_condition_1d = thermal_boundary_condition_data<float, 1u>;

using material_data_1d = material_data<data_mock, float, 1u>;
using material_data_2d = material_data<data_mock, float, 2u>;

using thermal_material_data_1d = thermal_material_data<float, 1u>;
using thermal_material_data_2d = thermal_material_data<float, 2u>;

}