#pragma once

#include <solvers/solver_2d/base/matrix_assembler_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I, class J>
class mass_matrix : public matrix_assembler_2d<T, I, J, 2> {
    using _base = matrix_assembler_2d<T, I, J, 2>;
    using block_t = metamath::types::square_matrix<T, 2>;

    static constexpr bool Symmetric = true;
    static constexpr size_t DoF = 2zu;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t>& theories,
                                const std::vector<bool>& is_inner);

public:
    explicit mass_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~mass_matrix() noexcept override = default;

    void compute(const std::vector<bool>& is_inner);
};

template<class T, class I, class J>
mass_matrix<T, I, J>::mass_matrix(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class J>
T mass_matrix<T, I, J>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * _base::mesh().jacobian(e, q);
    return integral;
}

template<class T, class I, class J>
void mass_matrix<T, I, J>::create_matrix_portrait(const std::unordered_map<std::string, theory_t>& theories,
                                                              const std::vector<bool>& is_inner) {
    const size_t rows = DoF * _base::mesh().process_nodes().size();
    const size_t cols = DoF * _base::mesh().container().nodes_count();
    _base::matrix().inner().resize(rows, cols);
    _base::matrix().bound().resize(rows, cols);
    _base::init_shifts(theories, is_inner, Symmetric);
    _base::init_indices(theories, is_inner, Symmetric);
}

template<class T, class I, class J>
void mass_matrix<T, I, J>::compute(const std::vector<bool>& is_inner) {
    logger::info() << "Mass matrix assembly started" << std::endl;
    const std::unordered_map<std::string, theory_t> theories = local_theories(_base::mesh().container());
    create_matrix_portrait(theories, is_inner);
    _base::calc_coeffs(theories, is_inner, Symmetric,
        [this](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const T integral = integrate_basic_pair(e, i, j);
            return block_t{integral, 0, 0, integral};
        },
        [](const std::string&, const size_t, const size_t, const size_t, const size_t) constexpr noexcept { return block_t{}; }
    );
    logger::info() << "Mass matrix assembly started" << std::endl;
}

}