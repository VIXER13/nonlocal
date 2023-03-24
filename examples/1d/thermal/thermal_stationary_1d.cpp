#include "make_element.hpp"

#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "mesh_1d_utils.hpp"
#include "thermal_config_data.hpp"

namespace {

template<class T, template<class, size_t> class Physics>
std::shared_ptr<nonlocal::mesh::mesh_1d<T>> make_mesh(
    const std::vector<nonlocal::config::segment_data<T, Physics>>& materials,
    const size_t element_order, const size_t quadrature_order) {
    std::vector<nonlocal::mesh::segment_data<T>> segments(materials.size());
    std::vector<T> search_radii(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, segments.size()}) {
        segments[i] = nonlocal::mesh::segment_data<T>{
            .length = materials[i].length,
            .elements = materials[i].elements_count
        };
        search_radii[i] = materials[i].model.search_radius;
    }
    auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
        nonlocal::make_element<T>(nonlocal::element_1d_order_t(1), 
                                  nonlocal::quadrature_1d_order_t(1)), segments);
    mesh->find_neighbours(search_radii);
    return mesh;
}

template<class T>
nonlocal::thermal::parameters_1d<T> make_thermal_parameters(
    const std::vector<nonlocal::config::segment_data<T, nonlocal::config::thermal_material_data>>& materials) {
    nonlocal::thermal::parameters_1d<T> parameters(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        parameters[i] = {
            .model = {
                .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{materials[i].model.nonlocal_radius},
                .local_weight = materials[i].model.local_weight
            },
            .physical = {
                .conductivity = materials[i].physical.conductivity,
                .capacity = materials[i].physical.capacity,
                .density = materials[i].physical.density
            },
        };
    return parameters;
}

template<std::floating_point T>
std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>> make_boundary_condition(
    const nonlocal::config::thermal_boundary_condition_data<T>& condition) {
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
    if (argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);

        using T = double;
        using I = int64_t;
        const nonlocal::config::stationary_thermal_1d_data<T> config_data{
            nonlocal::config::read_json(std::filesystem::path{argv[1]})
        };

        const auto mesh = make_mesh(config_data.materials, config_data.element_order, config_data.quadrature_order);
        auto solution = nonlocal::thermal::stationary_heat_equation_solver_1d<T, I>(
            mesh, make_thermal_parameters(config_data.materials),
            nonlocal::thermal::thermal_boundaries_conditions_1d<T>{
                make_boundary_condition<T>(config_data.boundaries.conditions.at("left")),
                make_boundary_condition<T>(config_data.boundaries.conditions.at("right"))
            },
            [value = config_data.equation.right_part](const T x) constexpr noexcept { return value; },
            config_data.equation.energy
        );
        
        std::cout << "integral = " << nonlocal::mesh::utils::integrate(*mesh, solution.temperature()) << std::endl;
        if (!std::filesystem::exists(config_data.save.folder()))
            std::filesystem::create_directories(config_data.save.folder());
        if (config_data.save.contains("temperature"))
            nonlocal::mesh::utils::save_as_csv(*mesh, solution.temperature(), config_data.save.path("temperature", ".csv"));
        if (config_data.save.contains("flux")) {
            solution.calc_flux();
            nonlocal::mesh::utils::save_as_csv(*mesh, solution.flux(), config_data.save.path("flux", ".csv"));
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}