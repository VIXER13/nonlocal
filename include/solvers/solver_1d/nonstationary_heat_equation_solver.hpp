#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "heat_capacity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_third_kind_1d.hpp"
#include "parameters_1d.hpp"

namespace nonlocal::thermal {

class _nonstationary_calc final {
    explicit _nonstationary_calc() noexcept = default;

    template<class T, class I>
    static void prepare_nonstationary_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                                            std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                            Eigen::SparseMatrix<T, Eigen::RowMajor, I>& C_inner,
                                            const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                            const T tau) {
        K_inner *= tau;
        K_inner += C_inner;
        if (boundary_condition.front().first == boundary_condition_t::FIRST_KIND)
            K_inner.coeffRef(0, 0) = 1;
        if (boundary_condition.back().first == boundary_condition_t::FIRST_KIND)
            K_inner.coeffRef(K_inner.rows()-1, K_inner.cols()-1) = 1;
        for(std::unordered_map<size_t, T>& matrix_part : K_bound)
            for(auto& [key, val] : matrix_part)
                val *= tau;
    }

    template<class T, class I, class Initi_Dist, class Right_Part>
    static void nonstationary_calc(const nonstationary_solver_parameters_1d<T>& solver_param,
                                   const Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                   const std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                   const Eigen::SparseMatrix<T, Eigen::RowMajor>& C_inner,
                                   const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                   const Init_Dist& init_dist, const Right_Part& right_part) {
       const T tau = (solver_param.time_interval.back() - solver_param.time_interval.front()) / solver_param.steps;
        Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh()->nodes_count());
        Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev(mesh()->nodes_count());
        Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr(mesh()->nodes_count());
        for(size_t i = 0; i < mesh()->nodes_count(); ++i)
            temperature_prev[i] = init_dist(mesh()->node_coord(i));
        Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
        if(sol_parameters.save_freq != std::numeric_limits<uintmax_t>::max())
            nonstationary_solver_logger(temperature_prev, sol_parameters, 0);
        for(size_t step = 1; step < sol_parameters.steps + 1; ++step) {
            f.setZero();
            const T t = sol_parameters.time_interval.front() + step * tau;
            const stationary_boundary_t stationary_bound = _base::nonstationary_boundary_to_stationary(bound_cond, t);
            _base::boundary_condition_second_kind(f, stationary_bound);
            _base::template integrate_right_part(f, [&right_part, t](const T x) { return right_part(t, x); });
            f *= tau;
            f += C_inner.template selfadjointView<Eigen::Upper>() * temperature_prev;
            _base::boundary_condition_first_kind(f, stationary_bound, K_bound);
            temperature_curr = solver.template solveWithGuess(f, temperature_prev);
            temperature_prev.swap(temperature_curr);
            if(step % sol_parameters.save_freq == 0)
                nonstationary_solver_logger(temperature_prev, sol_parameters, step);
        }
    }

public:
    template<class T, class I, class Init_Dist, class Right_Part, class Influence_Function>
    friend void nonstationary_heat_equation_solver_1d(const nonlocal_parameters_1d<T>& nonloc_param,
                                                      const heat_equation_parameters_1d<T>& equation_param,
                                                      const nonstationary_solver_parameters_1d<T>& solver_param,
                                                      const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                      const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                                      const Init_Dist& init_dist,
                                                      const Right_Part& right_part,
                                                      const Influence_Function& influence_function);
};

template<class T, class I, class Init_Dist, class Right_Part, class Influence_Function>
void nonstationary_heat_equation_solver_1d(const nonlocal_parameters_1d<T>& nonloc_param,
                                           const heat_equation_parameters_1d<T>& equation_param,
                                           const nonstationary_solver_parameters_1d<T>& solver_param,
                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                           const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                           const Init_Dist& init_dist,
                                           const Right_Part& right_part,
                                           const Influence_Function& influence_function) {
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(equation_param.lambda, {boundary_condition.front().first, boundary_condition.back().first},
                                      nonloc_param.p1, influence_function);

    heat_capacity_matrix_1d<T, I> capacity{mesh};
    capacity.calc_matrix(equation_param.c, equation_param.rho, {boundary_condition.front().first, boundary_condition.back().first});

    const T tau = (solver_param.time_interval.back() - solver_param.time_interval.front()) / solver_param.steps;
    prepare_nonstationary_matrix(conductivity.matrix_inner(), conductivity.matrix_bound(), capacity.matrix_inner(), boundary_condition, tau);

}

}

#endif