#ifndef STRUCTURAL_SOLVER_HPP
#define STRUCTURAL_SOLVER_HPP

#include <functional>
#include <algorithm>
#include <omp.h>
#include "finite_element_solver_base.hpp"
#include "../../Eigen/Eigen/Dense"
#include "../../Eigen/Eigen/Sparse"
//#include "Eigen/PardisoSupport"

namespace nonlocal::structural {

enum component : size_t { _11, _22, _12 };

enum class boundary_t : uint8_t {
    DISPLACEMENT = uint8_t(boundary_type::FIRST_KIND),
    PRESSURE     = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
using bound_cond = boundary_condition<T, boundary_t, 2>;

template<class Type>
struct parameters {
    Type nu = 0, // Коэффициент Пуассона
         E  = 0; // Модуль Юнга
};

template<class T>
struct distributed_load {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    std::array<std::function<T(const std::array<T, 2>&)>, 2> 
        func = { [](const std::array<T, 2>&) noexcept { return 0; },
                 [](const std::array<T, 2>&) noexcept { return 0; } };
};

template<class T, class I>
class structural_solver : protected finite_element_solver_base<T, I> {
    using _base = finite_element_solver_base<T, I>;
    using typename _base::Finite_Element_2D_Ptr;
    using typename _base::component;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::quad_shift;
    using _base::quad_coord;
    using _base::jacobian;
    using _base::first_node;
    using _base::last_node;

    std::array<T, 3> _D;
    //std::vector<std::vector<I>> _nodes_neighbors;

    // Матрица Гука, которая имеет следующий портрет:
    // arr[0] arr[1]   0
    // arr[1] arr[0]   0
    //   0      0    arr[2]
    static std::array<T, 3> hooke_matrix(const parameters<T>& params) noexcept {
        return {             params.E / (1 - params.nu*params.nu),
                 params.nu * params.E / (1 - params.nu*params.nu),
                       0.5 * params.E / (1 + params.nu) };
    }

    template<bool Proj, bool Approx>
    T integrate_loc(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) const {
        T integral = 0;
        static constexpr size_t k = Proj ^ Approx;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) / jacobian(quad_shift) *
                        (_D[k] * _base::template dNd< Proj>(e, i, q, quad_shift) * _base::template dNd< Approx>(e, j, q, quad_shift) +
                         _D[2] * _base::template dNd<!Proj>(e, i, q, quad_shift) * _base::template dNd<!Approx>(e, j, q, quad_shift));
        return integral;
    }

    template<bool Proj, bool Approx, class Influence_Function>
    T integrate_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,  size_t shiftL,
                       const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL,
                       const Influence_Function& influence_function) const {
        T integral = 0;
        const size_t sub_shift = shiftNL;
        static constexpr size_t k = Proj ^ Approx;
        for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL) {
            T inner_int_x = 0, inner_int_y = 0;
            for(size_t qNL = 0, shiftNL = sub_shift; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL) {
                const T influence_weight = eNL->weight(qNL) * influence_function(quad_coord(shiftL), quad_coord(shiftNL));
                inner_int_x += influence_weight * _base::template dNd< Approx>(eNL, jNL, qNL, shiftNL);
                inner_int_y += influence_weight * _base::template dNd<!Approx>(eNL, jNL, qNL, shiftNL);
            }
            integral += eL->weight(qL) * (_D[k] * inner_int_x * _base::template dNd< Proj>(eL, iL, qL, shiftL) +
                                          _D[2] * inner_int_y * _base::template dNd<!Proj>(eL, iL, qL, shiftL));
        }
        return integral;
    }

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K, Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                const T p1, const std::vector<bool>& inner_nodes) const {
        std::vector<std::unordered_set<I>> inner_portrait(K.rows()),
                                           bound_portrait(K_bound.rows());

        const auto indexator = [&inner_nodes, &inner_portrait, &bound_portrait, shift = 2 * first_node()](const I row, const I col) {
            if (inner_nodes[row] && inner_nodes[col]) {
                if (row <= col)
                    inner_portrait[row - shift].insert(col);
            } else if (row != col) {
                if (!inner_nodes[col])
                    bound_portrait[row - shift].insert(col);
            } else
                inner_portrait[row - shift].insert(col);
        };

        const auto structural_indexator = [&indexator](const I row, const I col) {
            indexator(row + I(X), col + I(X));
            indexator(row + I(X), col + I(Y));
            indexator(row + I(Y), col + I(X));
            indexator(row + I(Y), col + I(Y));
        };

        if (p1 > _base::MAX_LOCAL_WEIGHT) {
            _base::template mesh_run_loc(
                [this, &structural_indexator] (const size_t el, const size_t i, const size_t j) {
                    structural_indexator(2 * mesh().node_number(el, i), 2 * mesh().node_number(el, j));
                });
        } else {
            _base::template mesh_run_nonloc(
                [this, &structural_indexator](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
                    structural_indexator(2 * mesh().node_number(elL, iL), 2 * mesh().node_number(elNL, jNL));
                });
        }

        _base::convert_portrait(K, inner_portrait);
        _base::convert_portrait(K_bound, bound_portrait);
    }

    template<class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K, Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                     const T p1, const Influence_Function& influence_fun, const std::vector<bool>& inner_nodes) const {
        const auto filler_loc =
            [this, &K, &K_bound, &inner_nodes, p1, shift = 2 * first_node()]
            (const size_t el, const size_t i, const size_t j, const component proj, const component approx, const auto& integrate_rule) {
                const I row = 2 * mesh().node_number(el, i) + I(proj),
                        col = 2 * mesh().node_number(el, j) + I(approx);
                if (inner_nodes[row] && inner_nodes[col]) {
                    if (row <= col)
                        K.coeffRef(row - shift, col) += p1 * integrate_rule(mesh().element_2d(el), i, j, quad_shift(el));
                } else if (row != col) {
                    if (!inner_nodes[col])
                        K_bound.coeffRef(row - shift, col) += p1 * integrate_rule(mesh().element_2d(el), i, j, quad_shift(el));
                } else
                    K.coeffRef(row - shift, col) = 1;
            };

        _base::template mesh_run_loc(
            [this, &filler_loc](const size_t el, const size_t i, const size_t j) {
#define SIGNATURE const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift
                filler_loc(el, i, j, X, X, [this](SIGNATURE) { return integrate_loc<X, X>(e, i, j, quad_shift); });
                filler_loc(el, i, j, X, Y, [this](SIGNATURE) { return integrate_loc<X, Y>(e, i, j, quad_shift); });
                filler_loc(el, i, j, Y, X, [this](SIGNATURE) { return integrate_loc<Y, X>(e, i, j, quad_shift); });
                filler_loc(el, i, j, Y, Y, [this](SIGNATURE) { return integrate_loc<Y, Y>(e, i, j, quad_shift); });
#undef SIGNATURE
            });

        if (p1 < _base::MAX_LOCAL_WEIGHT) {
            const auto filler_nonloc =
                [this, &K, &K_bound, &inner_nodes, &influence_fun, p2 = 1 - p1, shift = 2 * first_node()]
                (const size_t elL, const size_t iL, const size_t elNL, const size_t jNL, const component proj, const component approx, const auto& integrate_rule) {
                    const I row = 2 * mesh().node_number(elL,  iL ) + I(proj),
                            col = 2 * mesh().node_number(elNL, jNL) + I(approx);
                    if (inner_nodes[row] && inner_nodes[col]) {
                        if (row <= col)
                            K.coeffRef(row - shift, col) += p2 * integrate_rule(mesh().element_2d(elL ), iL,  quad_shift(elL ),
                                                                                mesh().element_2d(elNL), jNL, quad_shift(elNL), influence_fun);
                    } else if (row != col)
                        if (!inner_nodes[col])
                            K_bound.coeffRef(row - shift, col) += p2 * integrate_rule(mesh().element_2d(elL ), iL,  quad_shift(elL ),
                                                                                      mesh().element_2d(elNL), jNL, quad_shift(elNL), influence_fun);
                };

            _base::template mesh_run_nonloc(
                [this, &filler_nonloc](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
#define SIGNATURE const Finite_Element_2D_Ptr& eL, const size_t iL, size_t shiftL, const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL, const Influence_Function& influence_function
                    filler_nonloc(elL, iL, elNL, jNL, X, X, [this](SIGNATURE) { return integrate_nonloc<X, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, X, Y, [this](SIGNATURE) { return integrate_nonloc<X, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, Y, X, [this](SIGNATURE) { return integrate_nonloc<Y, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, Y, Y, [this](SIGNATURE) { return integrate_nonloc<Y, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
#undef SIGNATURE
                }
            );
        }
    }

    template<class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K, Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                       const std::vector<bound_cond<T>>& bounds_cond, const T p1, const Influence_Function& influence_fun) {
        std::vector<bool> inner_nodes(2*mesh().nodes_count(), true);
        _base::template boundary_nodes_run([this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type(comp) == boundary_t::DISPLACEMENT)
                    inner_nodes[2 * mesh().node_number(b, el, i) + comp] = false;
        });

        double time = omp_get_wtime();
        create_matrix_portrait(K, K_bound, p1, inner_nodes);
        std::cout << "create_matrix_portrait: " << omp_get_wtime() - time << std::endl;

        time = omp_get_wtime();
        calc_matrix(K, K_bound, p1, influence_fun, inner_nodes);
        std::cout << "calc coeffs: " << omp_get_wtime() - time << std::endl;

//        std::cout << Eigen::MatrixXd{K} << std::endl << std::endl;
//        std::cout << Eigen::MatrixXd{K_bound} << std::endl << std::endl;
    }

    template<class Vector>
    std::array<std::vector<std::array<T, 3>>, 2> strains_and_stress_loc(const Vector& displacement) const {
        std::vector<std::array<T, 3>> strain(mesh().nodes_count(), std::array<T, 3>{}),
                                      stress(mesh().nodes_count(), std::array<T, 3>{});
        for(size_t el = 0; el < mesh().elements_count(); ++el) {
            const auto& e = mesh().element_2d(el);
            for(size_t i = 0; i < e->nodes_count(); ++i) {
                std::array<T, 4> jacobi_matrix = {};
                for(size_t j = 0; j < e->nodes_count(); ++j) {
                    const std::array<T, 2>& node = mesh().node(mesh().node_number(el, j));
                    jacobi_matrix[0] += node[0] * e->Nxi (j, e->node(i));
                    jacobi_matrix[1] += node[0] * e->Neta(j, e->node(i));
                    jacobi_matrix[2] += node[1] * e->Nxi (j, e->node(i));
                    jacobi_matrix[3] += node[1] * e->Neta(j, e->node(i));
                }

                std::array<T, 3> strain_loc = {};
                for(size_t j = 0; j < e->nodes_count(); ++j) {
                    const size_t node = 2*mesh().node_number(el, j);
                    const T jac = jacobian(jacobi_matrix),
                            dx1 = ( jacobi_matrix[3] * e->Nxi(j, e->node(i)) - jacobi_matrix[2] * e->Neta(j, e->node(i))) / jac,
                            dx2 = (-jacobi_matrix[1] * e->Nxi(j, e->node(i)) + jacobi_matrix[0] * e->Neta(j, e->node(i))) / jac;
                    strain_loc[_11] += dx1 * displacement[node  ];
                    strain_loc[_22] += dx2 * displacement[node+1];
                    strain_loc[_12] += dx1 * displacement[node+1] + 
                                       dx2 * displacement[node  ];
                }

                stress[mesh().node_number(el, i)][_11] += _D[0] * strain_loc[_11] + _D[1] * strain_loc[_22];
                stress[mesh().node_number(el, i)][_22] += _D[1] * strain_loc[_11] + _D[0] * strain_loc[_22];
                stress[mesh().node_number(el, i)][_12] += _D[2] * strain_loc[_12];
                for(size_t j = 0; j < 3; ++j)
                    strain[mesh().node_number(el, i)][j] += strain_loc[j];
            }
        }

        for(size_t i = 0; i < mesh().nodes_count(); ++i) {
            const size_t repeating_count = _base::nodes_elements_map(i).size();
            strain[i][_11] /=     repeating_count;
            strain[i][_22] /=     repeating_count;
            strain[i][_12] /= 2 * repeating_count;
            stress[i][_11] /=     repeating_count;
            stress[i][_22] /=     repeating_count;
            stress[i][_12] /= 2 * repeating_count;
        }

        return {std::move(strain), std::move(stress)};
    }

    std::vector<std::array<T, 3>> approx_strain_in_quad(const std::vector<std::array<T, 3>>& strain) const {
        std::vector<std::array<T, 3>> strain_in_quad(quad_shift(mesh().elements_count()), std::array<T, 3>{});
#pragma omp parallel for default(none) shared(strain, strain_in_quad)
        for(size_t el = 0; el < mesh().elements_count(); ++el) {
            const auto& e = mesh().element_2d(el);
            for(size_t q = 0, shift = quad_shift(el); q < e->qnodes_count(); ++q, ++shift)
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 3; ++comp)
                        strain_in_quad[shift][comp] += strain[mesh().node_number(el, i)][comp] * e->qN(i, q);
        }
        return std::move(strain_in_quad);
    }

    template<class Influence_Function>
    void stress_nonloc(      std::vector<std::array<T, 3>>& stress, 
                       const std::vector<std::array<T, 3>>& strains,
                       const T p1, const Influence_Function& influence_fun) const {
        const T p2 = 1 - p1;
        const std::vector<std::array<T, 3>> strains_in_quad = approx_strain_in_quad(strains);
#pragma omp parallel for default(none) shared(stress, influence_fun)
        for(size_t node = 0; node < mesh().nodes_count(); ++node) {
            std::unordered_set<I> neighbors;
            for(const I elL : _base::nodes_elements_map(node))
                for(const I elNL : _base::neighbors(elL))
                    neighbors.insert(elNL);
            for(const I elNL : neighbors) {
                const auto& eNL = mesh().element_2d(mesh().element_2d_type(elNL));
                for(size_t q = 0, shift = quad_shift(elNL); q < eNL->qnodes_count(); ++q, ++shift) {
                    const T influence_weight = p2 * eNL->weight(q) * jacobian(shift) * influence_fun(quad_coord(shift), mesh().node(node));
                    stress[node][0] += influence_weight * (_D[0] * strains_in_quad[shift][0] + _D[1] * strains_in_quad[shift][1]);
                    stress[node][1] += influence_weight * (_D[1] * strains_in_quad[shift][0] + _D[0] * strains_in_quad[shift][1]);
                    stress[node][2] += influence_weight *  _D[2] * strains_in_quad[shift][2];
                }
            }
        }
    }

public:
    explicit structural_solver(const std::shared_ptr<mesh::mesh_info<T, I>>& mesh, const parameters<T>& params)
        : _base{mesh}
        , _D{hooke_matrix(params)} {}

    ~structural_solver() override = default;

    template<class Vector>
    void save_as_vtk(const std::string& path, const Vector& U,
                     const std::vector<std::array<T, 3>>& strain, const std::vector<std::array<T, 3>>& stress) const;

    template<class Influence_Function>
    Eigen::Matrix<T, Eigen::Dynamic, 1> stationary(
        const std::vector<bound_cond<T>> &bounds_cond, //const distributed_load<T>& right_part,
        const T p1, const Influence_Function& influence_fun);

    template<class Vector, class Influence_Function>
    std::array<std::vector<std::array<T, 3>>, 2> strains_and_stress(
        const Vector& displacement, const T p1, const Influence_Function& influence_fun) const;

    T calc_energy(const std::vector<std::array<T, 3>>& strain, 
                  const std::vector<std::array<T, 3>>& stress) const;
};

template<class T, class I>
template<class Vector>
void structural_solver<T, I>::save_as_vtk(const std::string& path, const Vector& U,
    const std::vector<std::array<T, 3>>& strain, const std::vector<std::array<T, 3>>& stress) const {
    static constexpr std::string_view data_type = std::is_same_v<T, float> ? "float" : "double";

    if(2 * mesh().nodes_count() != size_t(U.size()))
        throw std::domain_error{"2 * mesh().nodes_count() != U.size()."};

    std::ofstream fout{path};
    fout.precision(20);

    mesh().save_as_vtk(fout);

    fout << "POINT_DATA " << mesh().nodes_count() << '\n';
    fout << "VECTORS Displacement " << data_type << '\n';
    for(size_t i = 0; i < mesh().nodes_count(); ++i)
        fout << U[2*i] << ' ' << U[2*i+1] << " 0\n";

    static constexpr std::array<std::string_view, 3>
        strain_number = {"strain11", "strain22", "strain12"},
        stress_number = {"stress11", "stress22", "stress12"};

    for(size_t comp = 0; comp < 3; ++comp) {
        fout << "SCALARS " << strain_number[comp] << ' ' << data_type << " 1\n"
             << "LOOKUP_TABLE default\n";
        for(size_t i = 0; i < mesh().nodes_count(); ++i)
            fout << strain[i][comp] << '\n';
    }

    for(size_t comp = 0; comp < 3; ++comp) {
        fout << "SCALARS " << stress_number[comp] << ' ' << data_type << " 1\n"
             << "LOOKUP_TABLE default\n";
        for(size_t i = 0; i < mesh().nodes_count(); ++i)
            fout << stress[i][comp] << '\n';
    }

    fout << "SCALARS mises " << data_type << " 1\n"
         << "LOOKUP_TABLE default\n";
    for(size_t i = 0; i < mesh().nodes_count(); ++i)
        fout << std::sqrt(stress[i][_11] * stress[i][_11] +
                          stress[i][_22] * stress[i][_22] -
                          stress[i][_11] * stress[i][_22] +
                      3 * stress[i][_12] * stress[i][_12]) << '\n';
}

template<class T, class I>
template<class Influence_Function>
Eigen::Matrix<T, Eigen::Dynamic, 1> structural_solver<T, I>::stationary(
    const std::vector<bound_cond<T>> &bounds_cond, //const distributed_load<T>& right_part,
    const T p1, const Influence_Function& influence_fun) {
    double time = omp_get_wtime();
    const size_t rows = 2 * (last_node() - first_node()),
                 cols = 2 * mesh().nodes_count();
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K      (rows, cols),
                                               K_bound(rows, cols);
    create_matrix(K, K_bound, bounds_cond, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
    //integrate_right_part(mesh, right_part, f);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    _base::PETSc_solver(f, K);
    //std::cout << "x = " << std::endl;
    //VecView(x, nullptr);

//    time = omp_get_wtime();
//    //Eigen::PardisoLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::Lower> solver{K};
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
//    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
//    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;
//
//    return std::move(u);
    return f;
}

template<class T, class I>
template<class Vector, class Influence_Function>
std::array<std::vector<std::array<T, 3>>, 2> structural_solver<T, I>::strains_and_stress(
    const Vector& displacement, const T p1, const Influence_Function& influence_fun) const {
    auto [strain, stress] = strains_and_stress_loc(displacement);
    if(p1 < _base::MAX_LOCAL_WEIGHT) { // Нелокальная задача
        for(size_t i = 0; i < mesh().nodes_count(); ++i)
            for(size_t j = 0; j < 3; ++j)
                stress[i][j] *= p1;
        stress_nonloc(stress, strain, p1, influence_fun);
    }
    return {std::move(strain), std::move(stress)};
}

template<class T, class I>
T structural_solver<T, I>::calc_energy(const std::vector<std::array<T, 3>>& strain, 
                                       const std::vector<std::array<T, 3>>& stress) const {
    T integral = 0;
    for(size_t el = 0; el < mesh().elements_count(); ++el) {
        const auto& e = mesh().element_2d(el);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t q = 0, shift = quad_shift(el); q < e->qnodes_count(); ++q, ++shift)
                integral += e->weight(q) * e->qN(i, q) * jacobian(shift) * 
                            (    strain[mesh().node_number(el, i)][_11] * stress[mesh().node_number(el, i)][_11] +
                                 strain[mesh().node_number(el, i)][_22] * stress[mesh().node_number(el, i)][_22] + 
                             2 * strain[mesh().node_number(el, i)][_12] * stress[mesh().node_number(el, i)][_12]);
    }
    return 0.5 * integral;
}

}

#endif