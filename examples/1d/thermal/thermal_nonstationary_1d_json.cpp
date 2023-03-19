#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "json_parser.hpp"
#include <iostream>

namespace {

using T = double;
using I = int64_t;

void save_step(nonlocal::thermal::heat_equation_solution_1d<T>&& solution, const json_data& data, const uintmax_t step) {
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.temperature(), data.not_template.path_to_save_temperature + std::to_string(step) + ".csv");
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.calc_flux(), data.not_template.path_to_save_flux + std::to_string(step) + ".csv");
}

void save_info(const std::vector<nonlocal::equation_parameters<1, T, nonlocal::thermal::parameters_1d>>& parameters,
               const std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2>& boundary_condition,
               const json_data& data) {
    std::ofstream info_file;
    const std::string way_to_file = data.not_template.path_to_save_info + ".txt";
    std::cout << "Info about calculation will be writen in " << way_to_file << std::endl;
    info_file.open(way_to_file);

    // Есть мысль передавать в функцию только название папки, куда сохранять результат,
    // а всю информацию о расчете писать в специальный файл <calculation_info.txt>.
    // Это позволит не захламлять название расчета, при этом сохраняя всю информацию о нем.
    size_t num = 0;
    for (const auto& current_material : data.materials) {
    info_file << "------------------№"<< ++num << "-----------------" << std::endl;
    info_file << "conductivity = " << std::to_string(current_material.conductivity) << std::endl;
    info_file << "density = "      << std::to_string(current_material.density)      << std::endl;
    info_file << "capacity = "     << std::to_string(current_material.capacity)     << std::endl;
    info_file << "local_weight = " << std::to_string(current_material.local_weight) << std::endl;
    info_file << "-------------------------------------" << std::endl;
    }

    info_file.close();
}

auto get_segment_data(const json_data& data) {
    std::vector<nonlocal::mesh::segment_data<T>> res;
    for (const auto& current_material : data.materials) 
        res.push_back({current_material.length, current_material.elements});

    return res;       
}

std::vector<T> get_radii(const json_data& data) {
    std::vector<T> res;
    for (const auto& current_material : data.materials) 
        res.push_back(current_material.nonlocal_radius);

    return res;  
}

auto get_boundary_conditions(const json_data& data) {
    std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2> res;
    if (data.not_template.left_boundary_kind == boundary_kind_type::TEMPERATURE)
        res[0] = std::make_unique<nonlocal::thermal::temperature_1d<T>>(data.left_boundary_value);
    else
        res[0] = std::make_unique<nonlocal::thermal::flux_1d<T>>(data.left_boundary_value);
    if (data.not_template.right_boundary_kind == boundary_kind_type::TEMPERATURE)
        res[1] = std::make_unique<nonlocal::thermal::temperature_1d<T>>(data.right_boundary_value);
    else
        res[1] = std::make_unique<nonlocal::thermal::flux_1d<T>>(data.right_boundary_value);

    return res;
}

}

int main(const int argc, const char *const *const argv) {
    try {
        json_data data;
        get_data_from_json(data, argv[1]);
        std::cout.precision(3);
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
            nonlocal::make_element<T>((nonlocal::element_type)data.not_template.element_order),
            get_segment_data(data)
        );
        const std::vector<T> radii = get_radii(data);
        using thermal_parameter = nonlocal::equation_parameters<1, T, nonlocal::thermal::parameters_1d>;
        std::vector<thermal_parameter> parameters;
        for (const auto& current_material : data.materials)
            parameters.push_back(thermal_parameter{
                .model = {
                    .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{current_material.nonlocal_radius},
                    .local_weight = current_material.local_weight
                },
                .physical = {
                    .conductivity = current_material.conductivity,
                    .capacity = current_material.capacity,
                    .density = current_material.conductivity
                }
            });
        mesh->find_neighbours(radii);

        const std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2> boundary_condition = get_boundary_conditions(data);

        nonlocal::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, data.time_step};
        solver.compute(
            parameters,
            boundary_condition,
            [&data](const T) constexpr noexcept { return data.right_part; }
        );

        uint64_t steps_count = data.not_template.steps_count + 1;
        save_info(parameters, boundary_condition, data);
        save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, data, 0);
        for(const uintmax_t step : std::ranges::iota_view{1u, steps_count}) {
            solver.calc_step(boundary_condition, [&data](const T x) constexpr noexcept { return data.right_part; });
            if (!(step % data.not_template.save_frequent)) 
                save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, data, step);
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