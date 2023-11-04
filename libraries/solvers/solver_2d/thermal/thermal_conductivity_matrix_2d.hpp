#ifndef NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP
#define NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include <string>

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
    T integrate_loc(const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter,  const std::vector<T>& solution,
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
                       const std::vector<T>& solution,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner, const bool is_symmetric, const bool is_neumann);

    void integral_condition(const bool is_symmetric);

public:
    explicit thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~thermal_conductivity_matrix_2d() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner,
                 const bool is_symmetric = true, const bool is_neumann = false,
                 const std::optional<std::vector<T>>& solution = std::nullopt);
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
        integrator(q, el.weight(q) / mesh::jacobian(_base::mesh().jacobi_matrix(e, q)),
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
    T integral = T{0};
    switch(const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC:
        integrate_loc(e, i, j, [this, &integral, &conductivity, e](const size_t q, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            integral += factor * conductivity[X][X](_base::mesh().quad_coord(e, q)) * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]);
        });
    break;

    case material_t::ORTHOTROPIC:
        integrate_loc(e, i, j, [this, &integral, &conductivity, e](const size_t q, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            const std::array<T, 2>& qcoord = _base::mesh().quad_coord(e, q);
            integral += factor * (conductivity[X][X](qcoord) * dNi[X] * dNj[X] + conductivity[Y][Y](qcoord) * dNi[Y] * dNj[Y]);
        });
    break;

    case material_t::ANISOTROPIC:
        integrate_loc(e, i, j, [this, &integral, &conductivity, e](const size_t q, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            T integral_term = T{0};
            const std::array<T, 2>& qcoord = _base::mesh().quad_coord(e, q);
            for(const size_t row : std::ranges::iota_view{0u, 2u})
                for(const size_t col : std::ranges::iota_view{0u, 2u})
                    integral_term += conductivity[row][col](qcoord) * dNi[row] * dNj[col];
            integral += factor * integral_term;
        });
    break;

    default:
        unknown_material(parameter.material);
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(
    const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter,  const std::vector<T>& solution, const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    switch(const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: 
        integrate_loc(e, i, j, [this, &integral, &conductivity, &solution, e](const size_t q, const T factor, const std::array<T, 2>& dNi, const std::array<T, 2>& dNj) {
            const size_t qshift = _base::mesh().quad_shift(e) + q;
            const std::array<T, 2>& qcoord = _base::mesh().quad_coord(qshift);
            integral += factor * conductivity[X][X](qcoord, solution[qshift]) * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]);
        });
    break;

    case material_t::ORTHOTROPIC: 
        return std::numeric_limits<T>::quiet_NaN();
    break;

    case material_t::ANISOTROPIC: 
        return std::numeric_limits<T>::quiet_NaN();
    break;

    default:
        unknown_material(parameter.material);
    }
    return integral;
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
    T integral = T{0};
    switch (const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        const auto inner_integrator = [&influence, &conductivity](const size_t, const T weightNL, const std::array<T, 2>& qcoordL, const std::array<T, 2>& qcoordNL) {
            return weightNL * conductivity[X][X](qcoordNL) * influence(qcoordL, qcoordNL);
        };
        integrate_nonloc(eL, eNL, iL, jNL, inner_integrator,
        [&integral](const T weightL, const std::array<T, 2>& dNi, const std::array<T, 2>& inner_integral) {
            integral += weightL * (dNi[X] * inner_integral[X] + dNi[Y] * inner_integral[Y]);
        }); 
    } break;

    default:
        unknown_material(parameter.material);
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(
    const parameter_2d<T, coefficients_t::SOLUTION_DEPENDENT>& parameter, const Influence_Function& influence,
    const std::vector<T>& solution,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    T integral = T{0};
    switch (const auto& conductivity = parameter.conductivity; parameter.material) {
    case material_t::ISOTROPIC: {
        const auto inner_integrator = [this, &influence, &conductivity, &solution, eNL](const size_t qNL, const T weightNL, const std::array<T, 2>& qcoordL, const std::array<T, 2>& qcoordNL) {
            const size_t qshift = _base::mesh().quad_shift(eNL) + qNL;
            return weightNL * conductivity[X][X](qcoordNL, solution[qshift]) * influence(qcoordL, qcoordNL);
        };
        integrate_nonloc(eL, eNL, iL, jNL, inner_integrator,
        [&integral](const T weightL, const std::array<T, 2>& dNi, const std::array<T, 2>& inner_integral) {
            integral += weightL * (dNi[X] * inner_integral[X] + dNi[Y] * inner_integral[Y]);
        });
    } break;

    default:
        unknown_material(parameter.material);
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::create_matrix_portrait(
    const std::unordered_map<std::string, theory_t> theories, const std::vector<bool>& is_inner, const bool is_symmetric, const bool is_neumann) {
    const size_t rows = _base::mesh().process_nodes().size() + (is_neumann && parallel_utils::is_last_process());
    const size_t cols = _base::mesh().container().nodes_count() + is_neumann;
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    if (is_neumann) {
        for(const size_t row : std::views::iota(0u, rows))
            _base::matrix_inner().outerIndexPtr()[row + 1] = 1;
        if (!is_symmetric && parallel_utils::is_last_process())
            _base::matrix_inner().outerIndexPtr()[rows] = cols - 1;
    }
    _base::init_shifts(theories, is_inner, is_symmetric);
    static constexpr bool SORT_INDICES = false;
    _base::init_indices(theories, is_inner, is_symmetric, SORT_INDICES);
    if (is_neumann) {
        for(const size_t row : std::ranges::iota_view{0u, rows}) {
            const size_t index = _base::matrix_inner().outerIndexPtr()[row + 1] - 1;
            _base::matrix_inner().innerIndexPtr()[index] = _base::mesh().container().nodes_count();
        }
        if (!is_symmetric && parallel_utils::is_last_process()) 
            for(const size_t col : std::ranges::iota_view{0u, cols}) {
                const size_t index = _base::matrix_inner().outerIndexPtr()[rows - 1] + col;
                _base::matrix_inner().innerIndexPtr()[index] = col;
            }
    }
    utils::sort_indices(_base::matrix_inner());
    utils::sort_indices(_base::matrix_bound());
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integral_condition(const bool is_symmetric) {
    const auto process_nodes = _base::mesh().process_nodes();
#pragma omp parallel for default(none) shared(process_nodes, is_symmetric)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node - process_nodes.front(), _base::mesh().container().nodes_count());
        for(const I e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        if (!is_symmetric && parallel_utils::is_last_process())
            _base::matrix_inner().coeffRef(_base::matrix_inner().rows() - 1, node) = val;
    }
    if (!is_symmetric && parallel_utils::MPI_size() > 1 && parallel_utils::is_last_process()) {
#pragma omp parallel for default(none) shared(process_nodes, is_symmetric)
        for(size_t node = 0; node < process_nodes.front(); ++node) {
            T& val = _base::matrix_inner().coeffRef(_base::matrix_inner().rows() - 1, node);
            for(const I e : _base::mesh().elements(node))
                val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        }
    }
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner, 
                                                                 const bool is_symmetric, const bool is_neumann,
                                                                 const std::optional<std::vector<T>>& solution) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_inner, is_symmetric, is_neumann);
    if (is_neumann)
        integral_condition(is_symmetric);
    _base::calc_coeffs(theories, is_inner, is_symmetric,
        [this, &parameters, &solution](const std::string& group, const size_t e, const size_t i, const size_t j) {
            using enum coefficients_t;
            const auto& [model, physic] = parameters.at(group);
            if (const auto* const parameter = parameter_cast<CONSTANTS>(physic.get()); parameter)
                return model.local_weight * integrate_loc(*parameter, e, i, j);
            if (const auto* const parameter = parameter_cast<SPACE_DEPENDENT>(physic.get()); parameter)
                return model.local_weight * integrate_loc(*parameter, e, i, j);
            if (const auto* const parameter = parameter_cast<SOLUTION_DEPENDENT>(physic.get()); parameter)
                return model.local_weight * integrate_loc(*parameter, *solution, e, i, j);
            return std::numeric_limits<T>::quiet_NaN();
        },
        [this, &parameters, &solution](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            using enum coefficients_t;
            const auto& [model, physic] = parameters.at(group);
            const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
            if (const auto* const parameter = parameter_cast<CONSTANTS>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(*parameter, model.influence, eL, eNL, iL, jNL);
            if (const auto* const parameter = parameter_cast<SPACE_DEPENDENT>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(*parameter, model.influence, eL, eNL, iL, jNL);
            if (const auto* const parameter = parameter_cast<SOLUTION_DEPENDENT>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(*parameter, model.influence, *solution, eL, eNL, iL, jNL);
            return std::numeric_limits<T>::quiet_NaN();
        });
}

}

#endif