#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include "finite_element_solver_base_1d.hpp"
#include <fstream>

namespace nonlocal::heat {

template<class T>
struct equation_parameters final {
    T lambda   = T{1}; // Коэффициент теплопроводности
    T rho      = T{1}; // Плотность
    T c        = T{1}; // Теплоёмкость
    T p1       = T{1}; // Весовой параметр модели
    T r        = T{0}; // Радиус нелокального влияния
    T integral = T{0}; // Значение интеграла для задачи Неймана
};

template<class T>
struct solver_parameters final {
    std::string save_path; // Путь куда сохранять данные
    std::array<T, 2> time_interval = {0, 1};
    uintmax_t steps = 100;
    uintmax_t save_freq = 1; // Частота сохранения
    bool save_csv    = true; // Сохранять .csv файлы в формате (x1, x2, T)
    bool calc_energy = true; // Вычислять энергия при сохранении, иногда полезно для контроля расчёта
};

inline constexpr boundary_condition_t TEMPERATURE = boundary_condition_t::FIRST_KIND;
inline constexpr boundary_condition_t FLOW = boundary_condition_t::FIRST_KIND;
inline constexpr boundary_condition_t IMPEDANCE = boundary_condition_t::THIRD_KIND;

template<class T>
class heat_equation_solver_1d final : public finite_element_solver_base_1d<T> {
    using _base = finite_element_solver_base_1d<T>;
    using typename _base::stationary_boundary_t;
    using typename _base::nonstatinary_boundary_t;
    using _base::mesh;

    enum class matrix : bool {THERMAL_CONDUCTIVITY, HEAT_CAPACITY};

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                const stationary_boundary_t& bound_cond,
                                const bool neumann_task, const bool nonlocal_task) const;

    void neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner) const;

    template<matrix Type, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                     std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                     const stationary_boundary_t& bound_cond,
                     const equation_parameters<T>& parameters,
                     const bool nonlocal_task, const Influence_Function& influence_function) const;

    void prepare_nonstationary_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                      std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                      Eigen::SparseMatrix<T, Eigen::RowMajor>& C_inner,
                                      std::array<std::unordered_map<size_t, T>, 2>& C_bound,
                                      const stationary_boundary_t& bound_cond,
                                      const equation_parameters<T>& parameters, const T tau) const;

    template<class Init_Dist, class Right_Part>
    void nonstationary_calc(const solver_parameters<T>& sol_parameters,
                            const Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                            const std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                            const Eigen::SparseMatrix<T, Eigen::RowMajor>& C_inner,
                            const nonstatinary_boundary_t& bound_cond,
                            const Init_Dist& init_dist, const Right_Part& right_part) const;

    void nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                     const solver_parameters<T>& sol_parameters, const uintmax_t step) const;

public:
    explicit heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~heat_equation_solver_1d() override = default;

    template<class Right_Part, class Influence_Function>
    std::vector<T> stationary(const equation_parameters<T>& parameters,
                              const stationary_boundary_t& bound_cond,
                              const Right_Part& right_part,
                              const Influence_Function& influence_function) const;

    template<class Init_Dist, class Right_Part, class Influence_Function>
    void nonstationary(const solver_parameters<T>& sol_parameters,
                       const equation_parameters<T>& parameters,
                       const nonstatinary_boundary_t& bound_cond,
                       const Init_Dist& init_dist,
                       const Right_Part& right_part,
                       const Influence_Function& influence_function) const;
};

template<class T>
heat_equation_solver_1d<T>::heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : finite_element_solver_base_1d<T>{mesh} {}

template<class T>
T heat_equation_solver_1d<T>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q);
    return integral * mesh()->jacobian();
}

template<class T>
T heat_equation_solver_1d<T>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->nodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q);
    return integral * mesh()->jacobian();
}

template<class T>
T heat_equation_solver_1d<T>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qNxi(i, q) * el->qNxi(j, q);
    return integral / mesh()->jacobian();
}

template<class T>
template<class Influence_Function>
T heat_equation_solver_1d<T>::integrate_nonloc(const size_t eL, const size_t eNL,
                                               const size_t iL, const size_t jNL,
                                               const Influence_Function& influence_function) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t qL = 0; qL < el->qnodes_count(); ++qL) {
        T inner_integral = T{0};
        const T qcoordL = mesh()->quad_coord(eL, qL);
        for(size_t qNL = 0; qNL < el->qnodes_count(); ++qNL) {
            const T qcoordNL = mesh()->quad_coord(eNL, qNL);
            inner_integral += el->weight(qNL) * influence_function(qcoordL, qcoordNL) * el->qNxi(jNL, qNL);
        }
        integral += el->weight(qL) * el->qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<class T>
void heat_equation_solver_1d<T>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                        const stationary_boundary_t& bound_cond,
                                                        const bool neumann_task, const bool nonlocal_task) const {
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.outerIndexPtr()[row+1] = 1;
    _base::create_matrix_portrait(K_inner, bound_cond, nonlocal_task);
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.innerIndexPtr()[K_inner.outerIndexPtr()[row+1]-1] = mesh()->nodes_count();
}

template<class T>
void heat_equation_solver_1d<T>::neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner) const {
#pragma omp parallel for default(none) shared(K_inner)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node) {
        T& val = K_inner.coeffRef(node, mesh()->nodes_count());
        for(const auto& [e, i] : mesh()->node_elements(node).arr)
            if(e != std::numeric_limits<size_t>::max())
                val += integrate_basic(e, i);
    }
}

template<class T>
template<typename heat_equation_solver_1d<T>::matrix Type, class Influence_Function>
void heat_equation_solver_1d<T>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                             std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                             const stationary_boundary_t& bound_cond,
                                             const equation_parameters<T>& parameters,
                                             const bool nonlocal_task, const Influence_Function& influence_function) const {
    if constexpr (Type == matrix::THERMAL_CONDUCTIVITY)
        _base::template calc_matrix(K_inner, K_bound, bound_cond, nonlocal_task, influence_function,
            [this, factor = parameters.lambda * parameters.p1](const size_t e, const size_t i, const size_t j) {
                return factor * integrate_loc(e, i, j);
            },
            [this, factor = parameters.lambda * (T{1} - parameters.p1)]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
                return factor * integrate_nonloc(eL, eNL, iL, jNL, influence_function);
            });
    else if constexpr (Type == matrix::HEAT_CAPACITY)
        _base::template calc_matrix(K_inner, K_bound, bound_cond, nonlocal_task, influence_function,
            [this](const size_t e, const size_t i, const size_t j) { return integrate_basic_pair(e, i, j); },
            [](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) { return 0; });
}

template<class T>
void heat_equation_solver_1d<T>::prepare_nonstationary_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                              std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                                              Eigen::SparseMatrix<T, Eigen::RowMajor>& C_inner,
                                                              std::array<std::unordered_map<size_t, T>, 2>& C_bound,
                                                              const stationary_boundary_t& bound_cond,
                                                              const equation_parameters<T>& parameters, const T tau) const {
    C_bound[0].clear();
    C_bound[1].clear();
    C_inner *= parameters.rho * parameters.c;
    K_inner *= tau;
    K_inner += C_inner;
    if (bound_cond[0].first == boundary_condition_t::FIRST_KIND)
        K_inner.coeffRef(0, 0) == 1;
    if (bound_cond[1].first == boundary_condition_t::FIRST_KIND)
        K_inner.coeffRef(K_inner.rows()-1, K_inner.cols()-1) = 1;
    for(size_t b = 0; b < K_bound.size(); ++b)
        for(auto& [key, val] : K_bound[b])
            val *= tau;
}

template<class T>
template<class Init_Dist, class Right_Part>
void heat_equation_solver_1d<T>::nonstationary_calc(const solver_parameters<T>& sol_parameters,
                                                    const Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                    const std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                                    const Eigen::SparseMatrix<T, Eigen::RowMajor>& C_inner,
                                                    const nonstatinary_boundary_t& bound_cond,
                                                    const Init_Dist& init_dist, const Right_Part& right_part) const {
    const T tau = (sol_parameters.time_interval.back() - sol_parameters.time_interval.front()) / sol_parameters.steps;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh()->nodes_count()),
                                        temperature_prev(mesh()->nodes_count()),
                                        temperature_curr(mesh()->nodes_count());
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

template<class T>
void heat_equation_solver_1d<T>::nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                                             const solver_parameters<T>& sol_parameters, const uintmax_t step) const {
    std::cout << "step = " << step << std::endl;
    if(sol_parameters.save_csv) {
        std::ofstream csv{sol_parameters.save_path + std::to_string(step) + ".csv"};
        csv.precision(std::numeric_limits<T>::max_digits10);
        const T h = (mesh()->section().back() - mesh()->section().front()) / (mesh()->nodes_count() - 1);
        for(size_t i = 0; i < mesh()->nodes_count(); ++i)
            csv << mesh()->section().front() + i * h << ',' << temperature[i] << '\n';
    }
    //if(sol_parameters.calc_energy)
    //    std::cout << "Energy = " << _base::mesh_proxy()->integrate_solution(temperature) << std::endl;
}

template<class T>
template<class Right_Part, class Influence_Function>
std::vector<T> heat_equation_solver_1d<T>::stationary(const equation_parameters<T>& parameters,
                                                      const stationary_boundary_t& bound_cond,
                                                      const Right_Part& right_part,
                                                      const Influence_Function& influence_function) const {
    const bool neumann_task = bound_cond[0].first == boundary_condition_t::SECOND_KIND &&
                              bound_cond[1].first == boundary_condition_t::SECOND_KIND;
    if (neumann_task && bound_cond[0].second + bound_cond[1].second > 1e-5)
        throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};

    const size_t size = mesh()->nodes_count() + neumann_task;
    Eigen::SparseMatrix<T, Eigen::RowMajor> K_inner(size, size);
    std::array<std::unordered_map<size_t, T>, 2> K_bound;

    const bool nonlocal_task = parameters.p1 < 0.999;
    create_matrix_portrait(K_inner, bound_cond, neumann_task, nonlocal_task);
    calc_matrix<matrix::THERMAL_CONDUCTIVITY>(K_inner, K_bound, bound_cond, parameters, nonlocal_task, influence_function);
    if (neumann_task)
        neumann_task_col_fill(K_inner);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(size);
    _base::template integrate_right_part(f, right_part);
    _base::boundary_condition_second_kind(f, bound_cond);
    _base::boundary_condition_first_kind(f, bound_cond, K_bound);
    if (neumann_task)
        f[f.size()-1] = parameters.integral;

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh()->nodes_count())};
}

template<class T>
template<class Init_Dist, class Right_Part, class Influence_Function>
void heat_equation_solver_1d<T>::nonstationary(const solver_parameters<T>& sol_parameters,
                                               const equation_parameters<T>& parameters,
                                               const nonstatinary_boundary_t& bound_cond,
                                               const Init_Dist& init_dist,
                                               const Right_Part& right_part,
                                               const Influence_Function& influence_function) const {
    static constexpr bool NOT_NEUMANN_TASK = false;
    static constexpr bool LOCAL = false;
    const stationary_boundary_t stationary_bound = _base::nonstationary_boundary_to_stationary(bound_cond, sol_parameters.time_interval.front());

    Eigen::SparseMatrix<T, Eigen::RowMajor> K_inner(mesh()->nodes_count(), mesh()->nodes_count());
    std::array<std::unordered_map<size_t, T>, 2> K_bound;
    const bool nonlocal_task = parameters.p1 < 0.999;
    create_matrix_portrait(K_inner, stationary_bound, NOT_NEUMANN_TASK, nonlocal_task);
    calc_matrix<matrix::THERMAL_CONDUCTIVITY>(K_inner, K_bound, stationary_bound, parameters, nonlocal_task, influence_function);

    Eigen::SparseMatrix<T, Eigen::RowMajor> C_inner(mesh()->nodes_count(), mesh()->nodes_count());
    std::array<std::unordered_map<size_t, T>, 2> C_bound;
    create_matrix_portrait(C_inner, stationary_bound, NOT_NEUMANN_TASK, LOCAL);
    calc_matrix<matrix::HEAT_CAPACITY>(C_inner, C_bound, stationary_bound, parameters, LOCAL, influence_function);

    const T tau = (sol_parameters.time_interval.back() - sol_parameters.time_interval.front()) / sol_parameters.steps;
    prepare_nonstationary_matrix(K_inner, K_bound, C_inner, C_bound, stationary_bound, parameters, tau);
    nonstationary_calc(sol_parameters, K_inner, K_bound, C_inner, bound_cond, init_dist, right_part);
}

}

#endif