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

protected:
    static hooke_parameters to_hooke(const parameters_2d<T>& parameters);

    static block_t calc_block(const hooke_matrix<T>& hooke, const block_t& integral) noexcept;
    static void add_to_integral(block_t& integral, const std::array<T, 2>& wdN, const std::array<T, 2>& dN) noexcept;
    block_t integrate_loc(const hooke_matrix<T>& hooke, const size_t e, const size_t i, const size_t j) const;
    // std::array<T, 4> integrate_nonloc(const hooke_parameter<T>& parameter,
    //                                   const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner);

public:
    explicit stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~stiffness_matrix() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner);
};

template<class T, class I, class J>
stiffness_matrix<T, I, J>::stiffness_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class J>
stiffness_matrix<T, I, J>::hooke_parameters stiffness_matrix<T, I, J>::to_hooke(const parameters_2d<T>& parameters) {
    hooke_parameters params;
    for(const auto& [group, equation_parameters] : parameters)
        params[group] = {.model = equation_parameters.model, .physical = equation_parameters.physical.hooke()};
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
    for(const size_t i : std::ranges::iota_view{0u, 2u})
        for(const size_t j : std::ranges::iota_view{0u, 2u})
            integral[i][j] += wdN[i] * dN[j];
}

template<class T, class I, class J>
stiffness_matrix<T, I, J>::block_t stiffness_matrix<T, I, J>::integrate_loc(const hooke_matrix<T>& hooke, const size_t e, const size_t i, const size_t j) const {
    block_t integral = {};
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
        using namespace metamath::functions;
        const T weight = el.weight(q) / mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
        add_to_integral(integral, weight * _base::mesh().derivatives(e, i, q), _base::mesh().derivatives(e, j, q));
    }
    return calc_block(hooke, integral);
}

// template<class T, class I, class J>
// template<class Influence_Function>
// std::array<T, 4> stiffness_matrix<T, I, J>::integrate_nonloc(const std::array<T, 3>& D,
//                                                                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
//                                                                         const Influence_Function& influence_function) const {
//     const auto& elL            = _base::mesh().element_2d(eL ),
//               & elNL           = _base::mesh().element_2d(eNL);
//           auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
//           auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
//     const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
//     const auto  dNdNL_begin    = _base::mesh_proxy()->dNdX(eNL, jNL);
//     std::array<T, 4> pairs = {};
//     for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
//         auto dNdNL    = dNdNL_begin;
//         auto qcoordNL = qcoordNL_begin;
//         std::array<T, 2> inner_int = {};
//         for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
//             const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
//             inner_int[X] += influence_weight * (*dNdNL)[X];
//             inner_int[Y] += influence_weight * (*dNdNL)[Y];
//         }
//         const std::array<T, 2> wdNdL = {elL->weight(qL) * (*dNdL)[X], elL->weight(qL) * (*dNdL)[Y]};
//         add_to_pair(pairs, wdNdL, inner_int);
//     }
//     return calc_block(D, pairs);
// }

// template<class T, class I, class J>
// template<class Influence_Function>
// void stiffness_matrix<T, I, J>::calc_matrix(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner) {
//     const theory_t theory = p1 < MAX_NONLOCAL_WEIGHT<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
//     const size_t rows = 2 * (_base::mesh_proxy()->last_node() - _base::mesh_proxy()->first_node()),
//                  cols = 2 * _base::mesh().nodes_count();
//     _base::matrix_inner().resize(rows, cols);
//     _base::matrix_bound().resize(rows, cols);
//     if (theory == theory_t::LOCAL)
//         _base::template create_matrix_portrait<theory_t::LOCAL>(is_inner);
//     else if (theory == theory_t::NONLOCAL)
//         _base::template create_matrix_portrait<theory_t::NONLOCAL>(is_inner);
//     _base::template calc_matrix(is_inner, theory, influence_fun,
//         [this, p1, &D](const size_t e, const size_t i, const size_t j) {
//             using namespace metamath::functions;
//             return p1 * integrate_loc(D, e, i, j);
//         },
//         [this, p2 = 1 - p1, &D](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
//             using namespace metamath::functions;
//             return p2 * integrate_nonloc(D, eL, eNL, iL, jNL, influence_function);
//         }
//     );
// }

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                                       const std::vector<bool>& is_inner) {
    const size_t rows = 2 * _base::mesh().process_nodes().size();
    const size_t cols = 2 * _base::mesh().container().nodes_count();
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    _base::init_shifts(theories, is_inner);
    _base::init_indices(theories, is_inner);
}

template<class T, class I, class J>
void stiffness_matrix<T, I, J>::compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_inner);
    const hooke_parameters hooke = to_hooke(parameters);
    _base::calc_coeffs(theories, is_inner,
        [this, &hooke](const std::string& group, const size_t e, const size_t i, const size_t j) {
            using namespace metamath::functions;
            const auto& parameter = hooke.at(group);
            return parameter.model.local_weight * integrate_loc(parameter.physical, e, i, j);
        },
        [this, &hooke](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& parameter = hooke.at(group);
            return nonlocal_weight(parameter.model.local_weight);
        }
    );
}

}

#endif