#pragma once

#include "mechanical_parameters_2d.hpp"

#include <solvers/solver_2d/base/matrix_assembler_2d.hpp>
#include <solvers/solver_2d/base/problem_settings.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I, class J>
class stiffness_matrix : public matrix_assembler_2d<T, I, J, 2> {
    using _base = matrix_assembler_2d<T, I, J, 2>;
    using hooke_parameter = equation_parameters<2, T, evaluated_hook_matrix_t>;
    using hooke_parameters = std::unordered_map<std::string, hooke_parameter>;
    using block_t = metamath::types::square_matrix<T, 2>;

    static constexpr bool NEUMANN = false;
    static constexpr bool SYMMETRIC = true;

protected:
    template<class Hooke>
    block_t integrate_local(const Hooke& hooke_matrix, const size_t e, const size_t i, const size_t j) const;
    template<class Hooke>
    block_t integrate_nonlocal(const Hooke& hooke_matrix, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>)>& influence,
                               const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const problem_settings& settings);

    T integrate_basic(const size_t e, const size_t i) const;
    void integral_condition();

public:
    explicit stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~stiffness_matrix() noexcept override = default;

    void compute(const evaluated_mechanical_parameters<T>& hooke, const problem_settings& settings, const assemble_part part = assemble_part::FULL);
};

template<class T, class I, class J>
stiffness_matrix<T, I, J>::stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class J>
template<class Hooke>
stiffness_matrix<T, I, J>::block_t stiffness_matrix<T, I, J>::integrate_local(const Hooke& hooke_matrix, const size_t e, const size_t i, const size_t j) const {
    block_t integral = {};
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
            const T shear1 = hooke[_66] * dNj[Y];
            const T shear2 = hooke[_66] * dNj[X];
            integral[X][X] += hooke[_11] * dNj[X] * wdNi[X] + shear1 * wdNi[Y];
            integral[X][Y] += hooke[_12] * dNj[Y] * wdNi[X] + shear2 * wdNi[Y];
            integral[Y][X] += hooke[_12] * dNj[X] * wdNi[Y] + shear1 * wdNi[X];
            integral[Y][Y] += std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ?
                              hooke[_11] * dNj[Y] * wdNi[Y] + shear2 * wdNi[X] :
                              hooke[_22] * dNj[Y] * wdNi[Y] + shear2 * wdNi[X];
        } else if constexpr (std::is_same_v<Hooke, evaluated_anisotropic_hook_matrix_t<T>>) {
            const T hooke26 = hooke[_26] * dNj[Y];
            const T hooke16 = hooke[_16] * dNj[X];
            const T shear1 = hooke16 + hooke[_66] * dNj[Y];
            const T shear2 = hooke26 + hooke[_66] * dNj[X];
            integral[X][X] += (hooke[_11] * dNj[X] + hooke[_16] * dNj[Y]) * wdNi[X] + shear1 * wdNi[Y];
            integral[X][Y] += (hooke[_12] * dNj[Y] + hooke16            ) * wdNi[X] + shear2 * wdNi[Y];
            integral[Y][X] += (hooke[_12] * dNj[X] + hooke26            ) * wdNi[Y] + shear1 * wdNi[X];
            integral[Y][Y] += (hooke[_22] * dNj[Y] + hooke[_26] * dNj[X]) * wdNi[Y] + shear2 * wdNi[X];
        } else
            static_assert(false, "Unsupported coefficient type.");
    }
    return integral;
}

template<class T, class I, class J>
template<class Hooke>
stiffness_matrix<T, I, J>::block_t stiffness_matrix<T, I, J>::integrate_nonlocal(
    const Hooke& hooke_matrix, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>)>& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    block_t integral = {};
    const size_t qshiftNL = _base::mesh().quad_shift(eNL);
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    for(const size_t qL : elL.qnodes()) {
        using namespace anisotropic_indices;
        using namespace metamath::operators;
        std::array<T, 6> inner_integral = {};
        const auto& qnodeL = _base::mesh().quad_coord(eL, qL);
        for(const size_t qNL : elNL.qnodes()) {
            const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshiftNL + qNL] : 
                                                       std::get<Constant>(hooke_matrix);
            const T weight = elNL.weight(qNL) * influence(qnodeL, _base::mesh().quad_coord(eNL, qNL));
            const auto wdNj = weight * _base::mesh().derivatives(eNL, jNL, qNL);
            if constexpr (std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ||
                          std::is_same_v<Hooke, evaluated_orthotropic_hook_matrix_t<T>>) {
                inner_integral[0] += hooke[_66] * wdNj[Y];
                inner_integral[1] += hooke[_66] * wdNj[X];
                inner_integral[2] += hooke[_11] * wdNj[X];
                inner_integral[3] += hooke[_12] * wdNj[Y];
                inner_integral[4] += hooke[_12] * wdNj[X];
                inner_integral[5] += std::is_same_v<Hooke, evaluated_isotropic_hook_matrix_t<T>> ?
                                     hooke[_11] * wdNj[Y] : hooke[_22] * wdNj[Y];
            } else if constexpr (std::is_same_v<Hooke, evaluated_anisotropic_hook_matrix_t<T>>) {
                const T hooke16 = hooke[_16] * wdNj[X];
                const T hooke26 = hooke[_26] * wdNj[Y];
                inner_integral[0] += hooke[_66] * wdNj[Y] + hooke16;
                inner_integral[1] += hooke[_66] * wdNj[X] + hooke26;
                inner_integral[2] += hooke[_11] * wdNj[X] + hooke[_16] * wdNj[Y];
                inner_integral[3] += hooke[_12] * wdNj[Y] + hooke16;
                inner_integral[4] += hooke[_12] * wdNj[X] + hooke26;
                inner_integral[5] += hooke[_22] * wdNj[Y] + hooke[_26] * wdNj[X];
            } else
                static_assert(false, "Unsupported coefficient type.");
        }
        const auto wdNi = elL.weight(qL) * _base::mesh().derivatives(eL, iL, qL);
        integral[X][X] += inner_integral[2] * wdNi[X] + inner_integral[0] * wdNi[Y];
        integral[X][Y] += inner_integral[3] * wdNi[X] + inner_integral[1] * wdNi[Y];
        integral[Y][X] += inner_integral[4] * wdNi[Y] + inner_integral[0] * wdNi[X];
        integral[Y][Y] += inner_integral[5] * wdNi[Y] + inner_integral[1] * wdNi[X];
    }
    return integral;
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                                       const problem_settings& settings) {
    const size_t cols = _base::cols() + NEUMANN;
    const size_t rows = _base::rows() == _base::cols() ?
                        cols : _base::rows() + (NEUMANN && parallel::is_last_process());
    _base::matrix().inner().resize(rows, cols);
    _base::matrix().bound().resize(rows, cols);
    if (NEUMANN) {
        for(const size_t row : std::views::iota(0u, _base::mesh().container().nodes_count()))
            _base::matrix().inner().outerIndexPtr()[2 * row + 1] = 1;
        _base::matrix().inner().outerIndexPtr()[2 * _base::mesh().container().nodes_count() + 1] = 1;
    }
    _base::init_shifts(theories, settings.is_inner_nodes, settings.is_symmetric());
    static constexpr bool Sort_Indices = false;
    _base::init_indices(theories, settings.is_inner_nodes, settings.is_symmetric(), Sort_Indices);
    if (NEUMANN) {
        for(const size_t row : std::ranges::iota_view{0u, _base::mesh().container().nodes_count()})
            _base::matrix().inner().innerIndexPtr()[_base::matrix().inner().outerIndexPtr()[2 * row + 1] - 1] = 2 * _base::mesh().container().nodes_count();
        _base::matrix().inner().innerIndexPtr()[_base::matrix().inner().nonZeros() - 1] = 2 * _base::mesh().container().nodes_count();
    }
    nonlocal::utils::sort_indices(_base::matrix().inner());
    nonlocal::utils::sort_indices(_base::matrix().bound());
    logger::info() << "Matrix portrait is formed" << std::endl;
}

template<class T, class I, class J>
T stiffness_matrix<T, I, J>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes())
        integral += el.weight(q) * el.qN(i, q) * _base::mesh().jacobian(e, q);
    return integral;
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::integral_condition() {
    const auto process_nodes = _base::mesh().process_nodes();
#pragma omp parallel for default(none) shared(process_nodes)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix().inner().coeffRef(2 * (node - process_nodes.front()), 2 * _base::mesh().container().nodes_count());
        for(const I e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
    }
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::compute(const evaluated_mechanical_parameters<T>& hooke, const problem_settings& settings, const assemble_part part) {
    logger::info() << "Stiffness matrix assembly started" << std::endl;
    const std::unordered_map<std::string, theory_t> theories = part == assemble_part::LOCAL ? 
                                                               local_theories(_base::mesh().container()) :
                                                               theories_types(hooke);
    create_matrix_portrait(theories, settings);
    if (NEUMANN)
        integral_condition();
    _base::calc_coeffs(theories, settings.is_inner_nodes, settings.is_symmetric(),
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