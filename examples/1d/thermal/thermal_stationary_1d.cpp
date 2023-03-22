#include "make_element.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "mesh_1d_utils.hpp"
#include "parse_thermal_1d.hpp"

namespace {

template<class T>
std::shared_ptr<nonlocal::mesh::mesh_1d<T>> make_mesh(const nonlocal::config::stationary_thermal_1d_data& config_data) {
    std::vector<nonlocal::mesh::segment_data<T>> segments(config_data.segments.size());
    std::vector<T> search_radii(config_data.segments.size());
    for(const size_t i : std::ranges::iota_view{0u, segments.size()}) {
        segments[i] = nonlocal::mesh::segment_data<T>{
            .length = config_data.segments[i].length,
            .elements = config_data.segments[i].elements_count
        };
        search_radii[i] = config_data.segments[i].model.search_radius;
    }
    auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
        nonlocal::make_element<T>(nonlocal::element_type(config_data.element_order)), segments);
    mesh->find_neighbours(search_radii);
    return mesh;
}

template<class T>
nonlocal::thermal::parameters_1d<T> make_parameters(const nonlocal::config::stationary_thermal_1d_data& config_data) {
    nonlocal::thermal::parameters_1d<T> parameters(config_data.segments.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        parameters[i] = {
            .model = {
                .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{config_data.segments[i].model.nonlocal_radius},
                .local_weight = config_data.segments[i].model.local_weight
            },
            .physical = {
                .conductivity = config_data.segments[i].physical.conductivity,
                .capacity = config_data.segments[i].physical.capacity,
                .density = config_data.segments[i].physical.density
            },
        };
    return parameters;
}

template<class T>
std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>> make_boundary_condition(
    const nonlocal::config::thermal_boundary_condition_1d& condition) {
    switch (condition.kind) {
    case nonlocal::thermal::boundary_condition_t::TEMPERATURE:
        return std::make_unique<nonlocal::thermal::temperature_1d<T>>(condition.temperature);

    case nonlocal::thermal::boundary_condition_t::FLUX:
        return std::make_unique<nonlocal::thermal::flux_1d<T>>(condition.flux);

    case nonlocal::thermal::boundary_condition_t::CONVECTION:
        return std::make_unique<nonlocal::thermal::convection_1d<T>>(condition.heat_transfer, condition.temperature);

    case nonlocal::thermal::boundary_condition_t::RADIATION:
        return std::make_unique<nonlocal::thermal::radiation_1d<T>>(condition.emissivity, T{0});

    case nonlocal::thermal::boundary_condition_t::COMBINED:
        return std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity, T{0});

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

}

int main(const int argc, const char *const *const argv) {
    if (argc < 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        const nonlocal::config::stationary_thermal_1d_data config_data{nonlocal::config::read_json_from_file(argv[1])};

        using T = double;
        using I = int64_t;
        const auto mesh = make_mesh<T>(config_data);
        auto solution = nonlocal::thermal::stationary_heat_equation_solver_1d<T, I>(
            mesh, make_parameters<T>(config_data),
            nonlocal::thermal::thermal_boundaries_conditions_1d<T>{
                make_boundary_condition<T>(config_data.equation.left),
                make_boundary_condition<T>(config_data.equation.right)
            },
            [value = config_data.equation.right_part](const T x) constexpr noexcept { return value; },
            config_data.equation.energy
        );
        
        std::cout << "integral = " << nonlocal::mesh::utils::integrate(*mesh, solution.temperature()) << std::endl;
        if (config_data.save.contains("temperature"))
          nonlocal::mesh::utils::save_as_csv(*mesh, solution.temperature(), config_data.save.path("temperature"));
        if (config_data.save.contains("flux"))
          nonlocal::mesh::utils::save_as_csv(*mesh, solution.temperature(), config_data.save.path("flux"));
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}