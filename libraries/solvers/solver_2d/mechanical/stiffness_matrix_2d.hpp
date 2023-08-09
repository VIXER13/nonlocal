#ifndef NONLOCAL_STIFFNESS_MATRIX_2D_HPP
#define NONLOCAL_STIFFNESS_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "mechanical_parameters_2d.hpp"

namespace nonlocal::mechanical {

template<class T, class I, class J>
class stiffness_matrix : public finite_element_matrix_2d<2, T, I, J> {
    using _base = finite_element_matrix_2d<2, T, I, J>;
    using hooke_parameter = equation_parameters<2, T, hooke_matrix>;
    using hooke_parameters = std::unordered_map<std::string, hooke_parameter>;
    using block_t = metamath::types::square_matrix<T, 2>;

    static constexpr bool SYMMETRIC = true;

protected:
    template<theory_t Theory>
    static hooke_parameters to_hooke(const parameters_2d<T>& parameters, const plane_t plane);

    static block_t calc_block(const hooke_matrix<T>& hooke, const block_t& integral) noexcept;
    static void add_to_integral(block_t& integral, const std::array<T, 2>& wdN, const std::array<T, 2>& dN) noexcept;
    block_t integrate_loc(const hooke_matrix<T>& hooke, const size_t e, const size_t i, const size_t j) const;
    block_t integrate_nonloc(const hooke_parameter& parameter, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner);

public:
    explicit stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~stiffness_matrix() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const plane_t plane, const std::vector<bool>& is_inner);
};

template<class T, class I, class J>
stiffness_matrix<T, I, J>::stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class J>
template<theory_t Theory>
stiffness_matrix<T, I, J>::hooke_parameters stiffness_matrix<T, I, J>::to_hooke(const parameters_2d<T>& parameters, const plane_t plane) {
    hooke_parameters params;
    for(const auto& [group, equation_parameters] : parameters) {
        const T factor = Theory == theory_t::LOCAL ?
                         equation_parameters.model.local_weight :
                         nonlocal_weight(equation_parameters.model.local_weight);
        using namespace metamath::functions;
        params[group] = {
            .model = equation_parameters.model,
            .physical = factor * equation_parameters.physical.hooke(plane)
        };
    }
    return params;
}

template<class T, class I, class J>
metamath::types::square_matrix<T, 2> stiffness_matrix<T, I, J>::calc_block(
    const hooke_matrix<T>& hooke, const metamath::types::square_matrix<T, 2>& integral) noexcept {
    return {
        hooke[0] * integral[X][X] + hooke[2] * integral[Y][Y],
        hooke[1] * integral[X][Y] + hooke[2] * integral[Y][X],
        hooke[1] * integral[Y][X] + hooke[2] * integral[X][Y],
        hooke[0] * integral[Y][Y] + hooke[2] * integral[X][X],
    };
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::add_to_integral(block_t& integral, const std::array<T, 2>& wdN, const std::array<T, 2>& dN) noexcept {
    for(const size_t i : std::ranges::iota_view{0u, wdN.size()})
        for(const size_t j : std::ranges::iota_view{0u, dN.size()})
            integral[i][j] += wdN[i] * dN[j];
}

template<class T, class I, class J>
stiffness_matrix<T, I, J>::block_t stiffness_matrix<T, I, J>::integrate_loc(
    const hooke_matrix<T>& hooke, const size_t e, const size_t i, const size_t j) const {
    block_t integral = {};
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
        using namespace metamath::functions;
        const T weight = el.weight(q) / mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
        add_to_integral(integral, weight * _base::mesh().derivatives(e, i, q), _base::mesh().derivatives(e, j, q));
    }
    return calc_block(hooke, integral);
}

template<class T, class I, class J>
stiffness_matrix<T, I, J>::block_t stiffness_matrix<T, I, J>::integrate_nonloc(
    const hooke_parameter& parameter, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    block_t integral = {};
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    for(const size_t qL : elL.qnodes()) {
        using namespace metamath::functions;
        std::array<T, 2> inner_integral = {};
        const auto influence = [&model = parameter.model, &qnodeL = _base::mesh().quad_coord(eL,  qL )](const std::array<T, 2>& qnodeNL) {
            return model.influence(qnodeL, qnodeNL);
        };
        for(const size_t qNL : elNL.qnodes())
            inner_integral += elNL.weight(qNL) * influence(_base::mesh().quad_coord(eNL, qNL)) * _base::mesh().derivatives(eNL, jNL, qNL);
        add_to_integral(integral, elL.weight(qL) * _base::mesh().derivatives(eL, iL, qL), inner_integral);
    }
    return calc_block(parameter.physical, integral);
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                                       const std::vector<bool>& is_inner) {
    const size_t rows = 2 * _base::mesh().process_nodes().size();
    const size_t cols = 2 * _base::mesh().container().nodes_count();
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    _base::init_shifts(theories, is_inner, SYMMETRIC);
    _base::init_indices(theories, is_inner, SYMMETRIC);
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::compute(const parameters_2d<T>& parameters, const plane_t plane, const std::vector<bool>& is_inner) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_inner);
    _base::calc_coeffs(theories, is_inner, SYMMETRIC,
        [this, hooke = to_hooke<theory_t::LOCAL>(parameters, plane)]
        (const std::string& group, const size_t e, const size_t i, const size_t j) {
            return integrate_loc(hooke.at(group).physical, e, i, j);
        },
        [this, hooke = to_hooke<theory_t::NONLOCAL>(parameters, plane)]
        (const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            return integrate_nonloc(hooke.at(group), eL, eNL, iL, jNL);
        }
    );
}

}

#endif