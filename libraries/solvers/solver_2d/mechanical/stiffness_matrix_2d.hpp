#pragma once

#include "mechanical_parameters_2d.hpp"

#include <solvers/solver_2d/base/matrix_assembler.hpp>
#include <solvers/solver_2d/base/problem_settings.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
class stiffness_matrix : public matrix_assembler_base<metamath::linear::square_matrix<T, 2>> {
    using _base = matrix_assembler_base<metamath::linear::square_matrix<T, 2>>;
    using hooke_parameter = equation_parameters<2, T, evaluated_hook_matrix_t>;
    using hooke_parameters = std::unordered_map<std::string, hooke_parameter>;

protected:
    template<class Hooke>
    metamath::linear::square_matrix<T, 2> integrate_local(
        const Hooke& hooke_matrix, const size_t e, const size_t i, const size_t j) const;
    template<class Hooke>
    metamath::linear::square_matrix<T, 2> integrate_nonlocal(
        const Hooke& hooke_matrix, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
        const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const problem_settings& settings);

public:
    explicit stiffness_matrix(const mesh::mesh_2d<T>& mesh);
    ~stiffness_matrix() noexcept override = default;

    void compute(const evaluated_mechanical_parameters<T>& hooke, const problem_settings& settings);
};

template<std::floating_point T>
stiffness_matrix<T>::stiffness_matrix(const mesh::mesh_2d<T>& mesh)
    : _base{mesh} {}

template<std::floating_point T>
template<class Hooke>
metamath::linear::square_matrix<T, 2> stiffness_matrix<T>::integrate_local(const Hooke& hooke_matrix, const size_t e, const size_t i, const size_t j) const {
    metamath::linear::square_matrix<T, 2> integral = {};
    const size_t qshift = _base::mesh().quad_shift(e);
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes()) {
        using namespace anisotropic_indices;
        using namespace metamath::operators;
        const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshift + q] : 
                                                   std::get<Constant>(hooke_matrix);
        const auto& dNj = _base::mesh().derivatives(e, j, q);
        const auto wdNi = (el.weight(q) / _base::mesh().jacobian(e, q)) * _base::mesh().derivatives(e, i, q);
        if constexpr (std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ||
                      std::is_same_v<Hooke, evaluated_orthotropic_hook_matrix_t<T>>) {
            const T shear1 = hooke[_66] * wdNi[Y];
            const T shear2 = hooke[_66] * wdNi[X];
            integral[X][X] += hooke[_11] * wdNi[X] * dNj[X] + shear1 * dNj[Y];
            integral[X][Y] += hooke[_12] * wdNi[X] * dNj[Y] + shear1 * dNj[X];
            integral[Y][X] += hooke[_12] * wdNi[Y] * dNj[X] + shear2 * dNj[Y];
            integral[Y][Y] += std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ?
                              hooke[_11] * wdNi[Y] * dNj[Y] + shear2 * dNj[X] :
                              hooke[_22] * wdNi[Y] * dNj[Y] + shear2 * dNj[X];
        } else if constexpr (std::is_same_v<Hooke, evaluated_anisotropic_hook_matrix_t<T>>) {
            const T hooke16 = hooke[_16] * wdNi[X];
            const T hooke26 = hooke[_26] * wdNi[Y];
            const T shear1 = hooke[_66] * wdNi[Y] + hooke16;
            const T shear2 = hooke[_66] * wdNi[X] + hooke26;
            integral[X][X] += (hooke[_11] * wdNi[X] + hooke[_16] * wdNi[Y]) * dNj[X] + shear1 * dNj[Y];
            integral[X][Y] += (hooke[_12] * wdNi[X] + hooke26             ) * dNj[Y] + shear1 * dNj[X];
            integral[Y][X] += (hooke[_12] * wdNi[Y] + hooke16             ) * dNj[X] + shear2 * dNj[Y];
            integral[Y][Y] += (hooke[_22] * wdNi[Y] + hooke[_26] * wdNi[X]) * dNj[Y] + shear2 * dNj[X];
        } else
            static_assert(false, "Unsupported coefficient type.");
    }
    return integral;
}

template<std::floating_point T>
template<class Hooke>
metamath::linear::square_matrix<T, 2> stiffness_matrix<T>::integrate_nonlocal(
    const Hooke& hooke_matrix, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    metamath::linear::square_matrix<T, 2> integral = {};
    const size_t qshiftNL = _base::mesh().quad_shift(eNL);
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    for(const size_t qL : elL.qnodes()) {
        using namespace anisotropic_indices;
        using namespace metamath::operators;
        const auto& qnodeL = _base::mesh().quad_coord(eL, qL);
        const auto wdNi = elL.weight(qL) * _base::mesh().derivatives(eL, iL, qL);
        for(const size_t qNL : elNL.qnodes()) {
            const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshiftNL + qNL] : 
                                                       std::get<Constant>(hooke_matrix);
            const T weight = elNL.weight(qNL) * influence(qnodeL, _base::mesh().quad_coord(eNL, qNL));
            const auto wdNj = weight * _base::mesh().derivatives(eNL, jNL, qNL);
            if constexpr (std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ||
                          std::is_same_v<Hooke, evaluated_orthotropic_hook_matrix_t<T>>) {
                const T shear1 = hooke[_66] * wdNi[Y];
                const T shear2 = hooke[_66] * wdNi[X];
                integral[X][X] += hooke[_11] * wdNi[X] * wdNj[X] + shear1 * wdNj[Y];
                integral[X][Y] += hooke[_12] * wdNi[X] * wdNj[Y] + shear1 * wdNj[X];
                integral[Y][X] += hooke[_12] * wdNi[Y] * wdNj[X] + shear2 * wdNj[Y];
                integral[Y][Y] += std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ?
                                  hooke[_11] * wdNi[Y] * wdNj[Y] + shear2 * wdNj[X] :
                                  hooke[_22] * wdNi[Y] * wdNj[Y] + shear2 * wdNj[X];
            } else if constexpr (std::is_same_v<Hooke, evaluated_anisotropic_hook_matrix_t<T>>) {
                const T hooke16 = hooke[_16] * wdNi[X];
                const T hooke26 = hooke[_26] * wdNi[Y];
                const T shear1 = hooke[_66] * wdNi[Y] + hooke16;
                const T shear2 = hooke[_66] * wdNi[X] + hooke26;
                integral[X][X] += (hooke[_11] * wdNi[X] + hooke[_16] * wdNi[Y]) * wdNj[X] + shear1 * wdNj[Y];
                integral[X][Y] += (hooke[_12] * wdNi[X] + hooke26             ) * wdNj[Y] + shear1 * wdNj[X];
                integral[Y][X] += (hooke[_12] * wdNi[Y] + hooke16             ) * wdNj[X] + shear2 * wdNj[Y];
                integral[Y][Y] += (hooke[_22] * wdNi[Y] + hooke[_26] * wdNi[X]) * wdNj[Y] + shear2 * wdNj[X];
            } else
                static_assert(false, "Unsupported coefficient type.");
        }
    }
    return integral;
}

template<std::floating_point T>
void stiffness_matrix<T>::create_matrix_portrait(const problem_settings& settings) {
    const size_t cols = _base::mesh().container().nodes_count();
    const size_t rows = _base::rows();
    _base::matrix().portrait.set_size(rows, cols);
    _base::init_shifts(settings);
    _base::init_indices(settings);
    logger::info() << "Matrix portrait is formed" << std::endl;
}

template<std::floating_point T>
void stiffness_matrix<T>::compute(const evaluated_mechanical_parameters<T>& hooke, const problem_settings& settings) {
    logger::info() << "Stiffness matrix assembly started" << std::endl;
    _base::matrix().clear();
    create_matrix_portrait(settings);
    _base::calc_coeffs(settings,
        [this, &hooke](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& [model, physic] = hooke.at(group);
            const auto integral = std::visit([this, e, i, j](const auto& hook) {
                return integrate_local(hook, e, i, j);
            }, physic.elastic);
            using namespace metamath::operators;
            return model.local_weight * integral;
        },
        [this, &hooke](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& [model, physic] = hooke.at(group);
            const auto integral = std::visit([this, &model, eL, eNL, iL, jNL](const auto& hook) {
                return integrate_nonlocal(hook, model.influence, eL, eNL, iL, jNL);
            }, physic.elastic);
            using namespace metamath::operators;
            return nonlocal::nonlocal_weight(model.local_weight) * integral;
        }
    );
    logger::info() << "Stiffness matrix assembly finished" << std::endl;
}

}