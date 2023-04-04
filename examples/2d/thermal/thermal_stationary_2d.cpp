#include "stationary_heat_equation_solver_2d.hpp"
#include "influence_functions_2d.hpp"
#include "mesh_container_2d_utils.hpp"

#include "thermal_config_data.hpp"

namespace nonlocal {

template<class T, class I>
auto make_mesh(const std::filesystem::path& path,
               const typename config::stationary_thermal_data<T, 2>::materials_t& materials) {
    auto mesh = std::make_shared<mesh::mesh_2d<T, I>>(path);
    std::unordered_map<std::string, T> radii;
    for(const auto& [name, material] : materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            radii.emplace(name, std::max(material.model.search_radius[0],
                                         material.model.search_radius[1]));
    mesh->find_neighbours(radii);
    return mesh;
}

template<class T>
thermal::parameters_2d<T> make_parameters(const typename config::stationary_thermal_data<T, 2>::materials_t& materials) {
    thermal::parameters_2d<T> parameters;
    for(const auto& [name, material] : materials) {
        auto& parameter = parameters[name] = equation_parameters<2, T, thermal::parameter_2d>{
            .model = {
                .influence = influence::polynomial_2d<T, 2, 1>{material.model.nonlocal_radius},
                .local_weight = material.model.local_weight
            },
            .physical = {
                .capacity = material.physical.capacity,
                .density = material.physical.density,
                .material = material.physical.material
            }
        };
        for(const size_t row : std::ranges::iota_view{0u, 2u})
            for(const size_t col : std::ranges::iota_view{0u, 2u})
                parameter.physical.conductivity[row][col] = material.physical.conductivity[2 * row + col];
    }
    return parameters;
}

template<std::floating_point T>
thermal::thermal_boundary_condition_2d<T> make_boundary_condition(
    const config::thermal_boundary_condition_data<T>& condition) {
    switch (condition.kind) {
    case thermal::boundary_condition_t::TEMPERATURE:
        return std::make_unique<thermal::temperature_2d<T>>(condition.temperature);

    case thermal::boundary_condition_t::FLUX:
        return std::make_unique<thermal::flux_2d<T>>(condition.flux);

    case thermal::boundary_condition_t::CONVECTION:
        return std::make_unique<thermal::convection_2d<T>>(condition.heat_transfer, condition.temperature);

    case thermal::boundary_condition_t::RADIATION:
        return std::make_unique<thermal::radiation_2d<T>>(condition.emissivity);

    case thermal::boundary_condition_t::COMBINED:
        return std::make_unique<thermal::combined_flux_2d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity);

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_2d<T> make_boundaries_conditions(
    const config::thermal_boundaries_conditions_data<T, 2>& conditions) {
    thermal::thermal_boundaries_conditions_2d<T> result;
    for(const auto& [name, condition] : conditions.conditions)
        result[name] = make_boundary_condition(condition);
    return result;
}

}

int main(const int argc, const char *const *const argv) {
#ifdef MPI_USE
    MPI_Init(&argc, &argv);
#endif

    if(argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    try {
        using T = double;
        using I = int64_t;
        const nonlocal::config::stationary_thermal_data<T, 2> config_data{nonlocal::config::read_json(std::filesystem::path{argv[1]})};
        std::cout.precision(config_data.other.get("precision", std::cout.precision()).asInt());

        const auto mesh = nonlocal::make_mesh<T, I>(config_data.mesh.mesh, config_data.materials);
        auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
            mesh, 
            nonlocal::make_parameters<T>(config_data.materials), 
            nonlocal::make_boundaries_conditions(config_data.boundaries), 
            [value = config_data.equation.right_part](const std::array<T, 2>& x) constexpr noexcept { return value; },
            config_data.equation.energy
        );
        solution.calc_flux();
        const auto& [TX, TY] = solution.flux();

        if (!std::filesystem::exists(config_data.save.folder()))
            std::filesystem::create_directories(config_data.save.folder());
        if (config_data.save.contains("config"))
            nonlocal::config::save_json(config_data.save.path("config", ".json"), config_data.to_json());
        if (config_data.save.contains("temperature"))
            nonlocal::mesh::utils::save_as_csv(config_data.save.path("temperature", ".csv"), mesh->container(), solution.temperature(), config_data.save.precision());
        if (config_data.save.contains("flux_x"))
            nonlocal::mesh::utils::save_as_csv(config_data.save.path("flux_x", ".csv"), mesh->container(), TX, config_data.save.precision());
        if (config_data.save.contains("flux_y"))
            nonlocal::mesh::utils::save_as_csv(config_data.save.path("flux_y", ".csv"), mesh->container(), TY, config_data.save.precision());
        if (config_data.save.contains("vtk"))
            solution.save_as_vtk(config_data.save.path("vtk", ".vtk"));
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}