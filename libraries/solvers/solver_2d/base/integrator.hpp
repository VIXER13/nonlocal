#ifndef NONLOCAL_INTEGRATOR_HPP
#define NONLOCAL_INTEGRATOR_HPP

#include "matrix_separator_base.hpp"

namespace nonlocal {

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
class integrator final : public matrix_separator_base<T, I>  {
    using _base = matrix_separator_base<T, I>;

    const mesh::mesh_container_2d<T, I>& _mesh;
    const Integrate_Loc& _integrate_loc;
    const Integrate_Nonloc& _integrate_nonloc;

public:
    explicit integrator(matrix_parts_t<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                        const std::vector<bool>& is_inner, const size_t node_shift,
                        const Integrate_Loc& integrate_loc, const Integrate_Nonloc& integrate_nonloc);
    ~integrator() noexcept override = default;

    void operator()(const std::string& group, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::integrator(
    matrix_parts_t<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
    const std::vector<bool>& is_inner, const size_t node_shift,
    const Integrate_Loc& integrate_loc, const Integrate_Nonloc& integrate_nonloc)
    : _base{matrix, is_inner, node_shift}
    , _mesh{mesh}
    , _integrate_loc{integrate_loc}
    , _integrate_nonloc{integrate_nonloc} {}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
void integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::operator()(
    const std::string& group, const size_t e, const size_t i, const size_t j) {
    const size_t row_glob = DoF * _mesh.node_number(e, i);
    const size_t col_glob = DoF * _mesh.node_number(e, j);
    std::optional<metamath::types::square_matrix<T, DoF>> block = std::nullopt;
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            if (const matrix_part part = _base::part(row, col); part != matrix_part::NO) {
                if (!block)
                    block = {_integrate_loc(group, e, i, j)};
                _base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col) += (*block)[row_loc][col_loc];
            }
        }
}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
void integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::operator()(
    const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    const size_t row_glob = DoF * _mesh.node_number(eL,  iL);
    const size_t col_glob = DoF * _mesh.node_number(eNL, jNL);
    std::optional<metamath::types::square_matrix<T, DoF>> block = std::nullopt;
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            if (const matrix_part part = _base::part(row, col); part != matrix_part::NO) {
                if (!block) {
                    block = {_integrate_nonloc(group, eL, eNL, iL, jNL)};
                    if (eL == eNL) {
                        const metamath::types::square_matrix<T, DoF> block_loc = {_integrate_loc(group, eL, iL, jNL)};
                        for(const size_t n : std::ranges::iota_view{0u, DoF})
                            for(const size_t m : std::ranges::iota_view{0u, DoF})
                                (*block)[n][m] += block_loc[n][m];
                    }
                }
                _base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col) += (*block)[row_loc][col_loc];
            }
        }
}

}

#endif