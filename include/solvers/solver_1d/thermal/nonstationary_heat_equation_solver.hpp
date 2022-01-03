#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_HPP

#include "thermal/thermal_conductivity_matrix_1d.hpp"
#include "heat_capacity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "convection_condition_1d.hpp"
#include "parameters_1d.hpp"
#include <iostream>
#include <fstream>

namespace nonlocal::thermal {

class _nonstationary_heat_equation_solver_1d final {
    explicit _nonstationary_heat_equation_solver_1d() noexcept = default;

    template<class T, class I>
    static void prepare_nonstationary_matrix(thermal_conductivity_matrix_1d<T, I>& K,
                                             const heat_capacity_matrix_1d<T, I>& C,
                                             const std::array<boundary_condition_t, 2> bound_cond,
                                             const T tau) {
        K.matrix_inner() *= tau;
        K.matrix_inner() += C.matrix_inner();
        if (bound_cond.front() == boundary_condition_t::TEMPERATURE)
            K.matrix_inner().coeffRef(0, 0) = T{1};
        if (bound_cond.back() == boundary_condition_t::TEMPERATURE)
            K.matrix_inner().coeffRef(K.matrix_inner().rows()-1, K.matrix_inner().cols()-1) = T{1};
        for(std::unordered_map<size_t, T>& matrix_part : K.matrix_bound())
            for(auto& [_, val] : matrix_part)
                val *= tau;
    }

    template<class T, class I, class Init_Dist, class Right_Part>
    static void nonstationary_calc(const nonstationary_solver_parameters_1d<T>& solver_param,
                                   const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                   const thermal_conductivity_matrix_1d<T, I>& K,
                                   const heat_capacity_matrix_1d<T, I>& C,
                                   const std::array<nonstatinary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                   const Init_Dist& init_dist, const Right_Part& right_part) {
       const T tau = (solver_param.time_interval.back() - solver_param.time_interval.front()) / solver_param.steps;
        Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count());
        Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev(mesh->nodes_count());
        Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr(mesh->nodes_count());
        for(const size_t i : std::views::iota(size_t{1}, mesh->nodes_count()))
            temperature_prev[i] = init_dist(mesh->node_coord(i));
        const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K.matrix_inner()};
        if(solver_param.save_freq != std::numeric_limits<uintmax_t>::max())
            nonstationary_solver_logger(temperature_prev, mesh, solver_param, 0);
        for(const uintmax_t step : std::views::iota(uintmax_t{1}, solver_param.steps + 1)) {
            f.setZero();
            const T t = solver_param.time_interval.front() + step * tau;
            const auto stationary_bound = to_stationary(boundary_condition, t);
            boundary_condition_second_kind_1d(f, stationary_bound, std::array{size_t{0}, size_t(f.size() - 1)});
            integrate_right_part(f, *mesh, [&right_part, t](const T x) { return right_part(t, x); });
            f *= tau;
            f += C.matrix_inner().template selfadjointView<Eigen::Upper>() * temperature_prev;
            boundary_condition_first_kind_1d(f, stationary_bound, K.matrix_bound());
            temperature_curr = solver.template solveWithGuess(f, temperature_prev);
            temperature_prev.swap(temperature_curr);
            if(step % solver_param.save_freq == 0)
                nonstationary_solver_logger(temperature_prev, mesh, solver_param, step);
        }
    }

    template<class T>
    static void nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                            const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                            const nonstationary_solver_parameters_1d<T>& solver_param,
                                            const uintmax_t step) {
        std::cout << "step = " << step << std::endl;
        if(solver_param.save_csv) {
            std::ofstream csv{solver_param.save_path.c_str() + std::to_string(step) + ".csv"};
            csv.precision(std::numeric_limits<T>::max_digits10);
            const T h = (mesh->section().back() - mesh->section().front()) / (mesh->nodes_count() - 1);
            for(const size_t i : std::views::iota(size_t{0}, mesh->nodes_count()))
                csv << mesh->section().front() + i * h << ',' << temperature[i] << '\n';
        }
    }

public:
    template<class T, class I, class Init_Dist, class Right_Part, class Influence_Function>
    friend void nonstationary_heat_equation_solver_1d(const heat_equation_parameters_1d<T>& equation_param,
                                                      const nonstationary_solver_parameters_1d<T>& solver_param,
                                                      const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                      const std::array<nonstatinary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                                      const Init_Dist& init_dist,
                                                      const Right_Part& right_part,
                                                      const T p1,
                                                      const Influence_Function& influence_function);
};

template<class T, class I, class Init_Dist, class Right_Part, class Influence_Function>
void nonstationary_heat_equation_solver_1d(const heat_equation_parameters_1d<T>& equation_param,
                                           const nonstationary_solver_parameters_1d<T>& solver_param,
                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                           const std::array<nonstatinary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                           const Init_Dist& init_dist,
                                           const Right_Part& right_part,
                                           const T p1,
                                           const Influence_Function& influence_function) {
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(equation_param.lambda, p1, influence_function, boundary_type(boundary_condition));
    convection_condition_1d(conductivity.matrix_inner(), boundary_type(boundary_condition), equation_param.alpha);

    heat_capacity_matrix_1d<T, I> capacity{mesh};
    capacity.calc_matrix(equation_param.c, equation_param.rho, boundary_type(boundary_condition));

    using _base = _nonstationary_heat_equation_solver_1d;
    const T tau = (solver_param.time_interval.back() - solver_param.time_interval.front()) / solver_param.steps;
    _base::prepare_nonstationary_matrix(conductivity, capacity, boundary_type(boundary_condition), tau);
    _base::nonstationary_calc(solver_param, mesh, conductivity, capacity,  boundary_condition, init_dist, right_part);
}

}

#endif