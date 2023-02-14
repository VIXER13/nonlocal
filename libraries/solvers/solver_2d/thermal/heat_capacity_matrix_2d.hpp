#ifndef NONLOCAL_HEAT_CAPACITY_MATRIX_2D_HPP
#define NONLOCAL_HEAT_CAPACITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class heat_capacity_matrix_2d : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t>& theories,
                                const std::vector<bool>& is_inner);

public:
    explicit heat_capacity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~heat_capacity_matrix_2d() noexcept override = default;

    void calc_matrix(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner);
};

template<class T, class I, class Matrix_Index>
heat_capacity_matrix_2d<T, I, Matrix_Index>::heat_capacity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class Matrix_Index>
T heat_capacity_matrix_2d<T, I, Matrix_Index>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.nodes_count()})
        integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
    return integral;
}

template<class T, class I, class Matrix_Index>
void heat_capacity_matrix_2d<T, I, Matrix_Index>::create_matrix_portrait(const std::unordered_map<std::string, theory_t>& theories,
                                                                         const std::vector<bool>& is_inner) {
    const size_t rows = _base::mesh().process_nodes().size();
    const size_t cols = _base::mesh().container().nodes_count();
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    _base::init_shifts(theories, is_inner);
    _base::init_indices(theories, is_inner);
}

template<class T, class I, class Matrix_Index>
void heat_capacity_matrix_2d<T, I, Matrix_Index>::calc_matrix(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner) {
    const auto theroires_setter = std::views::all(_base::mesh().container().groups_2d()) |
                                  std::views::transform([](const std::string& group) { return std::pair{group, theory_t::LOCAL}; });
    const std::unordered_map<std::string, theory_t> theories(theroires_setter.begin(), theroires_setter.end());
    create_matrix_portrait(theories, is_inner);
    _base::calc_coeffs(theories, is_inner,
        [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& parameter = parameters.at(group).physical;
            return parameter.density * parameter.capacity * integrate_basic_pair(e, i, j); 
        },
        [](const std::string&, const size_t, const size_t, const size_t, const size_t) constexpr noexcept { return T{0}; }
    );
}

}

#endif