#ifndef NONSTATIONARY_RELAX_TIME_HEAT_EQUATION_SOLVER_1D_HPP
#define NONSTATIONARY_RELAX_TIME_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "heat_capacity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

/*
    Class of the nonstationary 1d heat equation solver, which takes into account relaxation time.
    Problems:
        1)  Only works with the pulse heat on the left bound;
        2)  Param. _m (pulse heating intensity) can be changed only in the class;
        3)  Param. _relaxation_time can be changed only in the class
*/
namespace nonlocal::thermal {
template<class T, class I>
class nonstationary_relax_time_heat_equation_solver_1d final {
    heat_capacity_matrix_1d<T, I> _capacity;                        
    thermal_conductivity_matrix_1d<T, I> _conductivity;             
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> _left_matrix;       
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;                
    Eigen::Matrix<T, Eigen::Dynamic, 1> _integral_approx;           
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;          
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;          
    const int _m = 2;                                                            
    const T _time_step = T{1};                                      
    std::array<T, 2> _capacity_initial_values = {T(0), T(0)};       
    std::array<T, 2> _conductivity_initial_values = {T(0), T(0)};    
    std::array<T, 2> _left_matrix_initial_values = {T(0), T(0)};

    static std::array<T, 2> get_init_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
        return {matrix.coeff(0, 0), matrix.coeff(matrix.rows() - 1, matrix.cols() - 1)};
    }

    static void reset_to_init_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     const std::array<T, 2>& values) {
        matrix.coeffRef(0, 0) = values.front();
        matrix.coeffRef(matrix.rows() - 1, matrix.cols() - 1) = values.back();
    }
 
public:
    const T _relaxation_time = T{2};                                

    explicit nonstationary_relax_time_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step);

    const T time_step() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;

    template<class Init_Dist>
    void compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                 const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                 const Init_Dist& init_dist);

    template<class Right_Part>
    void calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                   const Right_Part& right_part, uint64_t time_iter);
};

template<class T, class I>
nonstationary_relax_time_heat_equation_solver_1d<T, I>::nonstationary_relax_time_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step)
    : _conductivity{mesh}
    , _capacity{mesh}
    , _left_matrix{_conductivity.matrix_inner()}   
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _integral_approx{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _time_step{time_step} {}

template<class T, class I>
const T nonstationary_relax_time_heat_equation_solver_1d<T, I>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_relax_time_heat_equation_solver_1d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I>
template<class Init_Dist>
void nonstationary_relax_time_heat_equation_solver_1d<T, I>::compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                                                          const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                          const Init_Dist& init_dist) {
    const std::array<bool, 2> is_first_kind = {
        bool(dynamic_cast<const temperature_1d<T>*>(boundaries_conditions.front().get())),
        bool(dynamic_cast<const temperature_1d<T>*>(boundaries_conditions.back ().get()))
    };
    _capacity.calc_matrix(parameters, is_first_kind);
    _conductivity.template calc_matrix(parameters, is_first_kind);
    convection_condition_1d(_conductivity.matrix_inner(), boundaries_conditions); 

    auto& capacity = _capacity.matrix_inner();
    if (_relaxation_time){
        capacity.coeffRef(0, 0) += capacity.coeffRef(0, 1);
        capacity.coeffRef(0, 1) = T{0};
        for (std::size_t i = 1; i < capacity.rows() - 1; i++) {
            capacity.coeffRef(i, i) += 2. * capacity.coeffRef(i, i + 1);
            capacity.coeffRef(i, i + 1) = T{0};
        }
        capacity.coeffRef(capacity.rows() - 1, capacity.cols() - 1) = capacity.coeffRef(0, 0);
    }                                            

    _left_matrix = _conductivity.matrix_inner();
    _left_matrix *= time_step() - _relaxation_time * (1. - exp(-time_step() / _relaxation_time));
    _left_matrix += _capacity.matrix_inner();
    reset_to_init_values(_left_matrix, {
        is_first_kind.front() ? T{1} : _left_matrix.coeffRef(0, 0),
        is_first_kind.back()  ? T{1} : _left_matrix.coeffRef(_left_matrix.rows() - 1, _left_matrix.cols() - 1)
    });
    _left_matrix_initial_values = get_init_values(_left_matrix);

    _conductivity_initial_values = get_init_values(_conductivity.matrix_inner());
    _capacity_initial_values = get_init_values(_capacity.matrix_inner());

    if constexpr (!std::is_same_v<Init_Dist, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        for(const size_t i : std::ranges::iota_view{0u, _conductivity.mesh().nodes_count()})
            _temperature_curr[i] = init_dist(_conductivity.mesh().node_coord(i));
}

template<class T, class I>
template<class Right_Part>
void nonstationary_relax_time_heat_equation_solver_1d<T, I>::calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                            const Right_Part& right_part, uint64_t time_iter) {
                                          
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    reset_to_init_values(_conductivity.matrix_inner(), _conductivity_initial_values);
    reset_to_init_values(_capacity.matrix_inner(), _capacity_initial_values);
    reset_to_init_values(_left_matrix, _left_matrix_initial_values);

    radiation_condition_1d(_capacity.matrix_inner(), boundaries_conditions, time_step());
    boundary_condition_second_kind_1d<T>(_right_part, boundaries_conditions);
    auto time = time_step() * time_iter;
    _right_part(0) += std::pow(_m, _m) / std::tgamma(_m) * std::pow(time, _m) * std::exp(-_m * time);  
    if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(_right_part, _conductivity.mesh(), right_part);
    _right_part *= time_step();
    _right_part += (_capacity.matrix_inner() - _relaxation_time * (1. - exp(-time_step() / _relaxation_time)) * _conductivity.matrix_inner()).template selfadjointView<Eigen::Upper>() * _temperature_prev;
    _right_part += _integral_approx;
    radiation_condition_1d(_conductivity.matrix_inner(), boundaries_conditions, time_step());
    boundary_condition_first_kind_1d(_right_part, _conductivity.matrix_bound(), boundaries_conditions);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_left_matrix};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);

    _integral_approx = std::exp(-time_step() / _relaxation_time) * ((_relaxation_time * (T{1} - std::exp(-time_step() / _relaxation_time)) * _conductivity.matrix_inner()).template selfadjointView<Eigen::Upper>() * (_temperature_curr - _temperature_prev) + _integral_approx);
}

}

#endif