#ifndef NONLOCAL_INTEGRATOR_HPP
#define NONLOCAL_INTEGRATOR_HPP

#include "matrix_separator_base.hpp"

namespace nonlocal {

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
class integrator final : public matrix_separator_base<T, I>  {
    using _base = matrix_separator_base<T, I>;
    const Integrate_Loc& _integrate_loc;
    const Integrate_Nonloc& _integrate_nonloc;

    bool check_block(const size_t row_block, const size_t col_block);
    void assemble(const metamath::types::square_matrix<T, DoF>& block, const size_t row_block, const size_t col_block);

public:
    explicit integrator(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift,
                        const Integrate_Loc& integrate_loc, const Integrate_Nonloc& integrate_nonloc);
    ~integrator() noexcept override = default;

    void operator()(const std::string& group, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::integrator(
    matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift,
    const Integrate_Loc& integrate_loc, const Integrate_Nonloc& integrate_nonloc)
    : _base{matrix, is_inner, node_shift}
    , _integrate_loc{integrate_loc}, _integrate_nonloc{integrate_nonloc} {}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
bool integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::check_block(const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF})
            if (_base::part(DoF * row_block + row_loc, DoF * col_block + col_loc) != matrix_part::NO)
                return true;
    return false;
}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
void integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::assemble(
    const metamath::types::square_matrix<T, DoF>& block, const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = DoF * row_block + row_loc;
            const size_t col = DoF * col_block + col_loc;
            _base::matrix(_base::part(row, col)).coeffRef(row - DoF *_base::node_shift, col) += block[row_loc][col_loc];
        }
}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
void integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::operator()(
    const std::string& group, const size_t e, const size_t i, const size_t j) {
    if (check_block(i, j))
        assemble({_integrate_loc(group, e, i, j)}, i, j);
}

template<size_t DoF, class T, class I, class Integrate_Loc, class Integrate_Nonloc>
void integrator<DoF, T, I, Integrate_Loc, Integrate_Nonloc>::operator()(
    const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    if (check_block(iL, jNL)) {
        metamath::types::square_matrix<T, DoF> block = {_integrate_nonloc(group, eL, eNL, iL, jNL)};
        if (eL == eNL) {
            const metamath::types::square_matrix<T, DoF> block_loc = {_integrate_loc(group, eL, iL, jNL)};
            for(const size_t row : std::ranges::iota_view{0u, DoF})
                for(const size_t col : std::ranges::iota_view{0u, DoF})
                    block[row][col] += block_loc[row][col];
        }
        assemble(block, iL, jNL);
    }
}

}

#endif