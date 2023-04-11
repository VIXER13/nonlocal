#ifndef NONLOCAL_THERMAL_CONFIG_DATA_HPP
#define NONLOCAL_THERMAL_CONFIG_DATA_HPP

#include "general_config_data.hpp"

#include <exception>

namespace nonlocal::config {

template<std::floating_point T>
struct thermal_equation_data final {
    T energy = 0;               // Used for Neumann problem
    T right_part = 0;
    T initial_distribution = 0; // Used for nonstationary and nonlinear problems

    explicit constexpr thermal_equation_data() noexcept = default;
    explicit thermal_equation_data(const Json::Value& equation);

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension = 0>
struct thermal_boundary_condition_data final {
    thermal::boundary_condition_t kind = thermal::boundary_condition_t::FLUX; // required
    T temperature = T{0};   // required if kind == TEMPERATURE or kind == CONVECTION
                            // used for TEMPERATURE condition, but if condition is CONVECTION used like ambient_temperature
    T flux = T{0};          // required if kind == FLUX
    T heat_transfer = T{0}; // required if kind == CONVECTION
    T emissivity = T{0};    // required if kind == RADIATION

    explicit constexpr thermal_boundary_condition_data() noexcept = default;
    explicit thermal_boundary_condition_data(const Json::Value& condition);

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension>
struct thermal_material_data;

template<std::floating_point T>
struct thermal_material_data<T, 1> final {
    T conductivity = T{1}; // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const Json::Value& physical);

    Json::Value to_json() const;
};

template<std::floating_point T>
class thermal_material_data<T, 2> final {
    void read_conductivity(const Json::Value& conduct);
    Json::Value save_conductivity() const;

public:
    material_t material = material_t::ISOTROPIC; // not json field
    std::array<T, 4> conductivity = {T{1}};      // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const Json::Value& physical);

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension>
using thermal_boundaries_conditions_data = boundaries_conditions_data<thermal_boundary_condition_data, T, Dimension>;



template<std::floating_point T, size_t Dimension>
struct stationary_thermal_data {
    using materials_t = std::conditional_t<
        Dimension == 1,
        std::vector<material_data<thermal_material_data, T, 1>>,
        std::unordered_map<std::string, material_data<thermal_material_data, T, Dimension>>
    >;

    Json::Value other;
    save_data save;
    mesh_data<Dimension> mesh;
    thermal_equation_data<T> equation;
    thermal_boundaries_conditions_data<T, Dimension> boundaries; // required
    materials_t materials;                                       // required

    explicit stationary_thermal_data(const Json::Value& value);
    virtual ~stationary_thermal_data() noexcept = default;

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension>
struct nonstationary_thermal_data final : public stationary_thermal_data<T, Dimension> {
    time_data<T> time;

    explicit nonstationary_thermal_data(const Json::Value& value);

    ~nonstationary_thermal_data() noexcept override = default;

    Json::Value to_json() const;
};

}

#include "thermal_config_data_impl.hpp"

#endif