#ifndef NONLOCAL_INTEGRATOR_HPP
#define NONLOCAL_INTEGRATOR_HPP

#include "matrix_separator_base.hpp"

namespace nonlocal {

template<class T, class I, size_t DoF, class Local_Integrator, class Nonlocal_Integrator>
class integrator final : public matrix_separator_base<T, I>  {
    using _base = matrix_separator_base<T, I>;
    using block_t = metamath::types::square_matrix<T, DoF>;

    const mesh::mesh_container_2d<T, I>& _mesh;
    const Local_Integrator& _local_integrator;
    const Nonlocal_Integrator& _nonlocal_integrator;

public:
    explicit integrator(finite_element_matrix<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                        const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric,
                        const Local_Integrator& local_integrator, const Nonlocal_Integrator& nonlocal_integrator);
    ~integrator() noexcept override = default;

    void operator()(const std::string& group, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<class T, class I, size_t DoF, class Local_Integrator, class Nonlocal_Integrator>
integrator<T, I, DoF, Local_Integrator, Nonlocal_Integrator>::integrator(
    finite_element_matrix<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
    const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric,
    const Local_Integrator& local_integrator, const Nonlocal_Integrator& nonlocal_integrator)
    : _base{matrix, is_inner, node_shift, is_symmetric}
    , _mesh{mesh}
    , _local_integrator{local_integrator}
    , _nonlocal_integrator{nonlocal_integrator} {}

template<class T, class I, size_t DoF, class Local_Integrator, class Nonlocal_Integrator>
void integrator<T, I, DoF, Local_Integrator, Nonlocal_Integrator>::operator()(
    const std::string& group, const size_t e, const size_t i, const size_t j) {
    const size_t row_glob = DoF * _mesh.node_number(e, i);
    const size_t col_glob = DoF * _mesh.node_number(e, j);
    std::optional<block_t> block = std::nullopt;
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF}) {
        T* val = nullptr;
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            if (const matrix_part part = _base::part(row, col); part != matrix_part::NO) {
                if (!block)
                    block = {_local_integrator(group, e, i, j)};
                if (part == matrix_part::BOUND)
                    _base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col) += (*block)[row_loc][col_loc];
                else {
                    val = val ?: &_base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col);
                    *val += (*block)[row_loc][col_loc];
                    ++val;
                }
            }
        }
    }
}

template<class T, class I, size_t DoF, class Local_Integrator, class Nonlocal_Integrator>
void integrator<T, I, DoF, Local_Integrator, Nonlocal_Integrator>::operator()(
    const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    const size_t row_glob = DoF * _mesh.node_number(eL,  iL);
    const size_t col_glob = DoF * _mesh.node_number(eNL, jNL);
    std::optional<block_t> block = std::nullopt;
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF}) {
        T* val = nullptr;
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            if (const matrix_part part = _base::part(row, col); part != matrix_part::NO) {
                if (!block) {
                    block = {_nonlocal_integrator(group, eL, eNL, iL, jNL)};
                    if (eL == eNL) {
                        using namespace metamath::functions;
                        (*block) += block_t{_local_integrator(group, eL, iL, jNL)};
                    }
                }
                if (part == matrix_part::BOUND)
                    _base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col) += (*block)[row_loc][col_loc];
                else {
                    val = val ?: &_base::matrix(part).coeffRef(row - DoF * _base::node_shift(), col);
                    *val += (*block)[row_loc][col_loc];
                    ++val;
                }
            }
        }
    }
}

}

#endif