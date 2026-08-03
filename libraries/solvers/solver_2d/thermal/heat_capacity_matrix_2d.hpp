#pragma once

#include <solvers/solver_2d/base/matrix_assembler.hpp>

namespace nonlocal::solver_2d::thermal {

template<class T>
class heat_capacity_matrix_2d : public matrix_assembler_base<T> {
    using _base = matrix_assembler_base<T>;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

    void create_matrix_portrait(const problem_settings& settings);

public:
    explicit heat_capacity_matrix_2d(const mesh::mesh_2d<T>& mesh);
    ~heat_capacity_matrix_2d() noexcept override = default;

    void calc_matrix(const parameters_2d<T>& parameters, const problem_settings& settings);
};

template<class T>
heat_capacity_matrix_2d<T>::heat_capacity_matrix_2d(const mesh::mesh_2d<T>& mesh)
    : _base{mesh} {}

template<class T>
T heat_capacity_matrix_2d<T>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * _base::mesh().jacobian(e, q);
    return integral;
}

template<class T>
void heat_capacity_matrix_2d<T>::create_matrix_portrait(const problem_settings& settings) {
    const size_t cols = _base::mesh().container().nodes_count();
    const size_t rows = _base::rows();
    _base::matrix().portrait.set_size(rows, cols);
    _base::init_shifts(settings);
    _base::init_indices(settings);
}

template<class T>
void heat_capacity_matrix_2d<T>::calc_matrix(const parameters_2d<T>& parameters, const problem_settings& settings) {
    const std::unordered_map<std::string, theory_t> theories = local_theories(_base::mesh().container());
    create_matrix_portrait(settings);
    _base::calc_coeffs(settings,
        [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& parameter = parameters.at(group).physical;
            return parameter.density * parameter.capacity * integrate_basic_pair(e, i, j); 
        },
        [](const std::string&, size_t, size_t, size_t, size_t) { return T{0}; }
    );
}

}