#ifndef NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP
#define NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "thermal_parameters_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class thermal_conductivity_matrix_2d : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;

    [[noreturn]] static void unknown_material(const material_t material);

protected:
    T integrate_basic(const size_t e, const size_t i) const;

    template<class Integrator>
    void integrate_loc(const size_t e, const size_t i, const size_t j, const Integrator& integrator) const;
    T integrate_loc(const parameter_2d<T, coefficients_t::CONSTANTS>& parameter, 
                    const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const parameter_2d<T, coefficients_t::SPACE_DEPENDENT>& parameter, 
                    const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter, 
                    const size_t e, const size_t i, const size_t j) const;
    
    template<class Integrator>
    std::array<T, 2> inner_integral(const size_t qnodes_count, const size_t eNL, const size_t jNL, const Integrator& integrator) const;
    template<class Inner_Integrator, class Integrator>
    void integrate_nonloc(const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                          const Inner_Integrator& inner_integrator, const Integrator& integrator) const;
    template<class Influence_Function>
    T integrate_nonloc(const parameter_2d<T, coefficients_t::CONSTANTS>& parameter, const Influence_Function& influence,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;
    template<class Influence_Function>
    T integrate_nonloc(const parameter_2d<T, coefficients_t::SPACE_DEPENDENT>& parameter, const Influence_Function& influence,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;
    template<class Influence_Function>
    T integrate_nonloc(const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter, const Influence_Function& influence,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner, const bool is_neumann);

    void neumann_problem_col_fill();

public:
    explicit thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~thermal_conductivity_matrix_2d() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner, const bool is_neumann = false);
};

template<class T, class I, class Matrix_Index>
thermal_conductivity_matrix_2d<T, I, Matrix_Index>::thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class Matrix_Index>
[[noreturn]] void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::unknown_material(const material_t material) {
    throw std::domain_error{"Unknown material type: " + std::to_string(std::underlying_type_t<material_t>(material))};
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes())
        integral += el.weight(q) * el.qN(i, q) * mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Integrator>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(const size_t e, const size_t i, const size_t j, const Integrator& integrator) const {
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes())
        integrator(q, el.weight(q) * mesh::jacobian(_base::mesh().jacobi_matrix(e, q)),
                   _base::mesh().derivatives(e, i, q), _base::mesh().derivatives(e, j, q));
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(
    const parameter_2d<T, coefficients_t::CONSTANTS>& parameter, const size_t e, const size_t i, const size_t j) const {
    switch(const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        T integral = T{0};
        integrate_loc(e, i, j, [&integral](const size_t, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            integral += factor * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]);
        });
        return conductivity[X][X] * integral;
    }

    case material_t::ORTHOTROPIC: {
        std::array<T, 2> integral_part = {};
        integrate_loc(e, i, j, [&integral_part](const size_t, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            integral_part[X] += factor * dNi[X] * dNj[X];
            integral_part[Y] += factor * dNi[Y] * dNj[Y];
        });
        return conductivity[X][X] * integral_part[X] + conductivity[Y][Y] * integral_part[Y];
    }

    case material_t::ANISOTROPIC: {
        metamath::types::square_matrix<T, 2> integral_part = {};
        integrate_loc(e, i, j, [&integral_part](const size_t, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            using namespace metamath::functions;
            const std::array<T, 2> fdNi = factor * dNi;
            for(const size_t row : std::ranges::iota_view{0u, 2u})
                for(const size_t col : std::ranges::iota_view{0u, 2u})
                    integral_part[row][col] += fdNi[row] * dNj[col];
        });
        return conductivity[X][X] * integral_part[X][X] + conductivity[X][Y] * integral_part[X][Y] +
               conductivity[Y][Y] * integral_part[Y][Y] + conductivity[Y][X] * integral_part[Y][X];
    }
    }
    unknown_material(parameter.material);
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(
    const parameter_2d<T, coefficients_t::SPACE_DEPENDENT>& parameter, const size_t e, const size_t i, const size_t j) const {
    switch(const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }

    case material_t::ORTHOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }

    case material_t::ANISOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }
    }
    unknown_material(parameter.material);
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(
    const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter, const size_t e, const size_t i, const size_t j) const {
    switch(const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }

    case material_t::ORTHOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }

    case material_t::ANISOTROPIC: {
        return std::numeric_limits<T>::quiet_NaN();
    }
    }
    unknown_material(parameter.material);
}

template<class T, class I, class Matrix_Index>
template<class Integrator>
std::array<T, 2> thermal_conductivity_matrix_2d<T, I, Matrix_Index>::inner_integral(
    const size_t qnodes_count, const size_t eNL, const size_t jNL, const Integrator& integrator) const {
    using namespace metamath::functions;
    std::array<T, 2> integral = {};
    for(const size_t qNL : std::ranges::iota_view{0u, qnodes_count})
        integral += integrator(qNL, _base::mesh().quad_coord(eNL, qNL)) * _base::mesh().derivatives(eNL, jNL, qNL);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Inner_Integrator, class Integrator>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
    const Inner_Integrator& inner_integrator, const Integrator& integrator) const {
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    for(const size_t qL : elL.qnodes()) {
        const auto callback = 
            [&inner_integrator, &elNL, &qcoordL = _base::mesh().quad_coord(eL, qL)](const size_t qNL, const std::array<T, 2>& qcoordNL) {
                return inner_integrator(qNL, elNL.weight(qNL), qcoordL, qcoordNL);
            };
        integrator(elL.weight(qL), _base::mesh().derivatives(eL, iL, qL), inner_integral(elNL.qnodes_count(), eNL, jNL, callback));
    }
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(
    const parameter_2d<T, coefficients_t::CONSTANTS>& parameter, const Influence_Function& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    const auto inner_integrator = [&influence](const size_t, const T weightNL, const std::array<T, 2>& qcoordL, const std::array<T, 2>& qcoordNL) {
        return weightNL * influence(qcoordL, qcoordNL);
    };
    switch (const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        T integral = T{0};
        integrate_nonloc(eL, eNL, iL, jNL, inner_integrator,
        [&integral](const T weightL, const std::array<T, 2>& dNi, const std::array<T, 2>& inner_integral) {
            integral += weightL * (dNi[X] * inner_integral[X] + dNi[Y] * inner_integral[Y]);
        });
        return conductivity[X][X] * integral;
    }
    
    case material_t::ORTHOTROPIC: {
        std::array<T, 2> integral_part = {};
        integrate_nonloc(eL, eNL, iL, jNL, inner_integrator,
        [&integral_part](const T weightL, const std::array<T, 2>& dNi, const std::array<T, 2>& inner_integral) {
            integral_part[X] += weightL * dNi[X] * inner_integral[X];
            integral_part[Y] += weightL * dNi[Y] * inner_integral[Y];
        });
        return conductivity[X][X] * integral_part[X] + conductivity[Y][Y] * integral_part[Y];
    }

    case material_t::ANISOTROPIC: {
        metamath::types::square_matrix<T, 2> integral_part = {};
        integrate_nonloc(eL, eNL, iL, jNL, inner_integrator,
        [&integral_part](const T weightL, const std::array<T, 2>& dNi, const std::array<T, 2>& inner_integral) {
            using namespace metamath::functions;
            const std::array<T, 2> wdNi = weightL * dNi;
            for(const size_t row : std::ranges::iota_view{0u, 2u})
                for(const size_t col : std::ranges::iota_view{0u, 2u})
                    integral_part[row][col] += wdNi[row] * inner_integral[col];
        });
        return conductivity[X][X] * integral_part[X][X] + conductivity[X][Y] * integral_part[X][Y] +
               conductivity[Y][Y] * integral_part[Y][Y] + conductivity[Y][X] * integral_part[Y][X];
    }
    }
    unknown_material(parameter.material);
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(
    const parameter_2d<T, coefficients_t::SPACE_DEPENDENT>& parameter, const Influence_Function& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    return 0;
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(
    const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter, const Influence_Function& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    return 0;
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                                                                const std::vector<bool>& is_inner, const bool is_neumann) {
    const size_t rows = _base::mesh().process_nodes().size() + (is_neumann && parallel_utils::MPI_rank() == parallel_utils::MPI_size() - 1);
    const size_t cols = _base::mesh().container().nodes_count() + is_neumann;
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    if (is_neumann)
        for(const size_t row : std::views::iota(0u, rows))
            _base::matrix_inner().outerIndexPtr()[row + 1] = 1;
    _base::init_shifts(theories, is_inner);
    static constexpr bool sort_indices = false;
    _base::init_indices(theories, is_inner, sort_indices);
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, rows})
            _base::matrix_inner().innerIndexPtr()[_base::matrix_inner().outerIndexPtr()[row + 1] - 1] = _base::mesh().container().nodes_count();
    utils::sort_indices(_base::matrix_inner());
    utils::sort_indices(_base::matrix_bound());
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::neumann_problem_col_fill() {
    const auto process_nodes = _base::mesh().process_nodes();
#pragma omp parallel for default(none) shared(process_nodes)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node - process_nodes.front(), _base::mesh().container().nodes_count());
        for(const I e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
    }
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::compute(const parameters_2d<T>& parameters,
                                                                 const std::vector<bool>& is_inner,
                                                                 const bool is_neumann) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_inner, is_neumann);
    _base::calc_coeffs(theories, is_inner,
        [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
            using enum coefficients_t;
            switch (const auto& [model, physic] = parameters.at(group); physic->type) {
            case CONSTANTS:
                return model.local_weight * integrate_loc(parameter_cast<CONSTANTS>(*physic), e, i, j);
            case SPACE_DEPENDENT:
                return model.local_weight * integrate_loc(parameter_cast<SPACE_DEPENDENT>(*physic), e, i, j);
            case SOLUTION_DEPENDENT:
                return model.local_weight * integrate_loc(parameter_cast<SOLUTION_DEPENDENT>(*physic), e, i, j);
            }
            return std::numeric_limits<T>::quiet_NaN();
        },
        [this, &parameters](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            using enum coefficients_t;
            switch (const auto& [model, physic] = parameters.at(group); physic->type) {
            case CONSTANTS:
                return nonlocal_weight(model.local_weight) * integrate_nonloc(parameter_cast<CONSTANTS>(*physic), model.influence, eL, eNL, iL, jNL);
            case SPACE_DEPENDENT:
                return nonlocal_weight(model.local_weight) * std::numeric_limits<T>::quiet_NaN();
            case SOLUTION_DEPENDENT:
                return nonlocal_weight(model.local_weight) * std::numeric_limits<T>::quiet_NaN();
            }
            return std::numeric_limits<T>::quiet_NaN();
        });
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif