#ifndef NONLOCAL_PARSE_THERMAL_1D_HPP
#define NONLOCAL_PARSE_THERMAL_1D_HPP

#include "config_utils.hpp"

#include "nonlocal_constants.hpp"

namespace nonlocal::config {

struct model_data final {
    double local_weight = 1.;
    double nonlocal_radius = 0.;
    double search_radius = 0.;

    explicit constexpr model_data() noexcept = default;
    explicit model_data(const Json::Value& model);
};

struct physical_data final {
    double conductivity = 1.;
    double capacity = 1.;
    double density = 1.;

    explicit constexpr physical_data() noexcept = default;
    explicit physical_data(const Json::Value& physical);
};

struct segment_data final {
    size_t elements_count = 100;
    double length = 1.;
    model_data model;
    physical_data physical;

    explicit constexpr segment_data() noexcept = default;
    explicit segment_data(const Json::Value& segment);
};

struct thermal_boundary_condition_1d final {
    thermal::boundary_condition_t kind = thermal::boundary_condition_t::FLUX;
    double temperature = 0.; // used for TEMPERATURE condition, but if condition is CONVECTION used like ambient_temperature
    double flux = 0.;
    double heat_transfer = 0.;
    double emissivity = 0.;

    static thermal::boundary_condition_t convert(const std::string& kind);

    explicit constexpr thermal_boundary_condition_1d() noexcept = default;
    explicit thermal_boundary_condition_1d(const Json::Value& condition);
};

struct thermal_equation_1d_data final {
    double energy = 0.; // Used for Neumann problem
    double right_part = 0.;
    thermal_boundary_condition_1d left;
    thermal_boundary_condition_1d right;

    explicit constexpr thermal_equation_1d_data() noexcept = default;
    explicit thermal_equation_1d_data(const Json::Value& equation);
};

struct stationary_thermal_1d_data {
    Json::Value other;
    save_data save;
    thermal_equation_1d_data equation;
    std::vector<segment_data> segments;
    size_t element_order = 1;

    static size_t convert(const std::string& order);

    explicit stationary_thermal_1d_data(const Json::Value& value);
    virtual ~stationary_thermal_1d_data() noexcept = default;
};

struct nonstationary_data final {
    double time_step = 0.01;
    double initial_time = 0.;
    double initial_distribution = 0.;
    uint64_t steps_cont = 100;
    uint64_t save_frequency = 1;

    explicit constexpr nonstationary_data() noexcept = default;
    explicit nonstationary_data(const Json::Value& nonstationary);
};

struct nonstationary_thermal_1d_data final : public stationary_thermal_1d_data {
    nonstationary_data nonstationary;

    explicit nonstationary_thermal_1d_data(const Json::Value& value);
    ~nonstationary_thermal_1d_data() noexcept override = default;
};

}

#endif