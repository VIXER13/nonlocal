#ifndef STRUCTURAL_SOLVER_HPP
#define STRUCTURAL_SOLVER_HPP

#include <tuple>
#include <functional>
#include <algorithm>
#include "omp.h"
#include "finite_element_solver_base.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/PardisoSupport"

namespace nonlocal::structural {

enum component : size_t { _11, _22, _12 };

enum class boundary_t : uint8_t {
    DISPLACEMENT = uint8_t(boundary_type::FIRST_KIND),
    PRESSURE     = uint8_t(boundary_type::SECOND_KIND)
};

template<class Type>
struct parameters {
    Type nu = 0, // Коэффициент Пуассона
         E  = 0; // Модуль Юнга
};

template<class T>
struct boundary_condition {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    std::array<std::function<T(const std::array<T, 2>&)>, 2> 
        func = { [](const std::array<T, 2>&) noexcept { return 0; },
                 [](const std::array<T, 2>&) noexcept { return 0; } };
    std::array<boundary_t, 2> type = { boundary_t::PRESSURE, boundary_t::PRESSURE };
};

template<class T>
struct distributed_load {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    std::array<std::function<T(const std::array<T, 2>&)>, 2> 
        func = { [](const std::array<T, 2>&) noexcept { return 0; },
                 [](const std::array<T, 2>&) noexcept { return 0; } };
};

template<class T, class I>
class structural_solver : protected finite_element_solver_base<T, I>
{
    using _base = finite_element_solver_base<T, I>;

    using typename _base::Finite_Element_1D_Ptr;
    using typename _base::Finite_Element_2D_Ptr;

    using _base::X;
    using _base::Y;
    using typename _base::component;
    using _base::MAX_LOCAL_WEIGHT;

    using _base::mesh;
    using _base::quad_shift;
    using _base::quad_coord;
    using _base::jacobi_matrix;

    using _base::jacobian;
    using _base::approx_quad_nodes_on_bound;
    using _base::approx_jacobi_matrices_on_bound;

    std::array<T, 3> _D;

    // Матрица Гука, которая имеет следующий портрет:
    // arr[0] arr[1]   0
    // arr[1] arr[0]   0
    //   0      0    arr[2]
    static std::array<T, 3> hooke_matrix(const parameters<T>& params) noexcept {
        return {             params.E / (1 - params.nu*params.nu),
                 params.nu * params.E / (1 - params.nu*params.nu),
                       0.5 * params.E / (1 + params.nu) };
    }

    template<bool Proj, bool Form>
    T integrate_loc(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) const {
        T integral = 0;
        static constexpr size_t k = Proj ^ Form;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) / jacobian(quad_shift) *
                        (_D[k] * _base::template dNd< Proj>(e, i, q, quad_shift) * _base::template dNd< Form>(e, j, q, quad_shift) +
                         _D[2] * _base::template dNd<!Proj>(e, i, q, quad_shift) * _base::template dNd<!Form>(e, j, q, quad_shift));
        return integral;
    }

    template<bool Proj, bool Form, class Influence_Function>
    T integrate_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,  size_t shiftL,
                       const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL,
                       const Influence_Function& influence_function) {
        T integral = 0;
        const size_t sub_shift = shiftNL;
        static constexpr size_t k = Proj ^ Form;
        for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL) {
            T inner_int_x = 0, inner_int_y = 0;
            for(size_t qNL = 0, shiftNL = sub_shift; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL) {
                const T influence_weight = eNL->weight(qNL) * influence_function(quad_coord(shiftL), quad_coord(shiftNL));
                inner_int_x += influence_weight * _base::template dNd< Form>(eNL, jNL, qNL, shiftNL);
                inner_int_y += influence_weight * _base::template dNd<!Form>(eNL, jNL, qNL, shiftNL);
            }
            integral += eL->weight(qL) * (_D[k] * inner_int_x * _base::template dNd< Proj>(eL, iL, qL, shiftL) +
                                          _D[2] * inner_int_y * _base::template dNd<!Proj>(eL, iL, qL, shiftL));
        }
        return integral;
    }

    std::array<std::vector<I>, 4> mesh_analysis(const std::vector<bool>& inner_nodes, const bool nonlocal) {
        std::vector<I> shifts_loc      (mesh().elements_count()+1, 0), 
                       shifts_bound_loc(mesh().elements_count()+1, 0),
                       shifts_nonloc, shifts_bound_nonloc;

        const auto counter_loc = 
            [this, &inner_nodes, &shifts_loc, &shifts_bound_loc]
            (const size_t el, const size_t i, const size_t j, const component proj, const component form) {
                const size_t row = 2 * mesh().node_number(el, i) + I(proj),
                             col = 2 * mesh().node_number(el, j) + I(form);
                if(row >= col) {
                    if(inner_nodes[row] && inner_nodes[col])
                        ++shifts_loc[el+1];
                    else if(row != col)
                        ++shifts_bound_loc[el+1];
                }
            };

        _base::template mesh_run_loc(
            [&counter_loc](const size_t el, const size_t i, const size_t j) {
                counter_loc(el, i, j, X, X);
                counter_loc(el, i, j, X, Y);
                counter_loc(el, i, j, Y, X);
                counter_loc(el, i, j, Y, Y);
            });

        shifts_loc[0] = std::count(inner_nodes.cbegin(), inner_nodes.cend(), false);
        for(size_t i = 1; i < shifts_loc.size(); ++i) {
            shifts_loc[i] += shifts_loc[i-1];
            shifts_bound_loc[i] += shifts_bound_loc[i-1];
        }

        if(nonlocal) {
            shifts_nonloc.resize(mesh().elements_count()+1, 0);
            shifts_bound_nonloc.resize(mesh().elements_count()+1, 0);

            const auto counter_nonloc =
                [this, &inner_nodes, &shifts_nonloc, &shifts_bound_nonloc]
                (const size_t elL, const size_t iL, const size_t elNL, const size_t jNL, const component proj, const component form) {
                    const size_t row = 2 * mesh().node_number(elL , iL ) + I(proj),
                                 col = 2 * mesh().node_number(elNL, jNL) + I(form);
                    if(row >= col) {
                        if(inner_nodes[row] && inner_nodes[col])
                            ++shifts_nonloc[elL+1];
                        else if(row != col)
                            ++shifts_bound_nonloc[elL+1];
                    }
                };

            _base::template mesh_run_nonloc(
                [&counter_nonloc](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
                    counter_nonloc(elL, iL, elNL, jNL, X, X);
                    counter_nonloc(elL, iL, elNL, jNL, X, Y);
                    counter_nonloc(elL, iL, elNL, jNL, Y, X);
                    counter_nonloc(elL, iL, elNL, jNL, Y, Y);
                });

            shifts_nonloc[0] = shifts_loc.back();
            shifts_bound_nonloc[0] = shifts_bound_loc.back();
            for(size_t i = 1; i < shifts_nonloc.size(); ++i) {
                shifts_nonloc[i] += shifts_nonloc[i-1];
                shifts_bound_nonloc[i] += shifts_bound_nonloc[i-1];
            }
        }

        return {std::move(shifts_loc), std::move(shifts_bound_loc), std::move(shifts_nonloc), std::move(shifts_bound_nonloc)};
    }

    template<class Influence_Function>
    std::array<std::vector<Eigen::Triplet<T, I>>, 2>
    triplets_fill(const std::vector<boundary_condition<T>> &bounds_cond, const T p1, const Influence_Function& influence_fun) {
        const bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
        std::vector<bool> inner_nodes(2*mesh().nodes_count(), true);
        _base::template boundary_nodes_run([this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT)
                    inner_nodes[2*mesh().node_number(b, el, i)+comp] = false;
        });

        auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(inner_nodes, nonlocal);

        const size_t triplets_count = nonlocal ? shifts_nonloc.back() + shifts_bound_nonloc.back()
                                               : shifts_loc.back()    + shifts_bound_loc.back();
        std::cout << "Triplets count: " << triplets_count << std::endl;
        std::vector<Eigen::Triplet<T, I>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                          triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
        for(size_t i = 0, j = 0; i < inner_nodes.size(); ++i)
            if(!inner_nodes[i])
                triplets[j++] = Eigen::Triplet<T, I>(i, i, 1.);

        const auto filler_loc =
            [this, &inner_nodes, &shifts_loc, &shifts_bound_loc, &triplets, &triplets_bound, p1]
            (const size_t el, const size_t i, const size_t j, const component proj, const component form, const auto& integrate_rule) {
                const I row = 2 * mesh().node_number(el, i) + I(proj),
                        col = 2 * mesh().node_number(el, j) + I(form);
                if(row >= col) {
                    T integral = p1 * integrate_rule(mesh().element_2d(el), i, j, quad_shift(el));
                    if(inner_nodes[row] && inner_nodes[col])
                        triplets[shifts_loc[el]++] = Eigen::Triplet<T, I>{row, col, integral};
                    else if(row != col)
                        triplets_bound[shifts_bound_loc[el]++] = inner_nodes[col] ? Eigen::Triplet<T, I>{col, row, integral} :
                                                                                    Eigen::Triplet<T, I>{row, col, integral};
                }
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

        if(nonlocal) {
            const auto filler_nonloc =
                [this, &inner_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc, &influence_fun, p2 = 1. - p1]
                (const size_t elL, const size_t iL, const size_t elNL, const size_t jNL, const component proj, const component form, const auto& integrate_rule) {
                    const I row = 2 * mesh().node_number(elL,  iL ) + I(proj),
                            col = 2 * mesh().node_number(elNL, jNL) + I(form);
                    if(row >= col) {
                        const T integral = p2 * integrate_rule(mesh().element_2d(elL ), iL,  quad_shift(elL ),
                                                               mesh().element_2d(elNL), jNL, quad_shift(elNL), influence_fun);
                        if(inner_nodes[row] && inner_nodes[col])
                            triplets[shifts_nonloc[elL]++] = Eigen::Triplet<T, I>{row, col, integral};
                        else if(row != col)
                            triplets_bound[shifts_bound_nonloc[elL]++] = inner_nodes[col] ? Eigen::Triplet<T, I>{col, row, integral} :
                                                                                            Eigen::Triplet<T, I>{row, col, integral};
                    }
                };

            _base::template mesh_run_nonloc(
                [this, &filler_nonloc](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
#define SIGNATURE const Finite_Element_2D_Ptr& eL,  const size_t iL,  size_t shiftL, const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL, const Influence_Function& influence_function
                    filler_nonloc(elL, iL, elNL, jNL, X, X, [this](SIGNATURE) { return integrate_nonloc<X, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, X, Y, [this](SIGNATURE) { return integrate_nonloc<X, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, Y, X, [this](SIGNATURE) { return integrate_nonloc<Y, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                    filler_nonloc(elL, iL, elNL, jNL, Y, Y, [this](SIGNATURE) { return integrate_nonloc<Y, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
#undef SIGNATURE
                });
        }

        return {std::move(triplets), std::move(triplets_bound)};
    }

    template<class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::ColMajor, I>& K, Eigen::SparseMatrix<T, Eigen::ColMajor, I>& K_bound,
                       const std::vector<boundary_condition<T>>& bounds_cond,const T p1, const Influence_Function& influence_fun) {
        const double time = omp_get_wtime();
        auto [triplets, triplets_bound] = triplets_fill(bounds_cond, p1, influence_fun);
        std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;

        K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
        triplets_bound.reserve(0);
        K.setFromTriplets(triplets.cbegin(), triplets.cend());
        std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
    }

    template<class Vector>
    void integrate_boundary_pressure(Vector& f, const std::vector<boundary_condition<T>>& bounds_cond) const {
        std::vector<std::array<T, 2>> quad_nodes, jacobi_matrices;
        for(size_t b = 0; b < bounds_cond.size(); ++b)
            if(bounds_cond[b].type[0] == boundary_t::PRESSURE || bounds_cond[b].type[1] == boundary_t::PRESSURE)
                for(size_t el = 0; el < mesh().elements_count(b); ++el) {
                    approx_quad_nodes_on_bound(quad_nodes, b, el);
                    approx_jacobi_matrices_on_bound(jacobi_matrices, b, el);
                    const auto& be = mesh().element_1d(b, el);
                    for(size_t i = 0; i < be->nodes_count(); ++i) 
                        for(size_t comp = 0; comp < 2; ++comp)
                            if(bounds_cond[b].type[comp] == boundary_t::PRESSURE)
                                f[2*mesh().node_number(b, el, i)+comp] += 
                                    _base::template integrate_boundary_gradient(be, i, quad_nodes, jacobi_matrices, bounds_cond[b].func[comp]);
                }
    }

    // Учёт граничных условий первого рода.
    template<class Vector>
    void displacement_on_boundary(Vector& f, const std::vector<boundary_condition<T>>& bounds_cond,
                                  const Eigen::SparseMatrix<T, Eigen::ColMajor, I>& K_bound) const {
        std::vector<std::vector<I>> kinematic_nodes(mesh().boundary_groups_count());
        _base::template boundary_nodes_run(
            [this, &bounds_cond, &kinematic_nodes](const size_t b, const size_t el, const size_t i) {
                if(bounds_cond[b].type[0] == boundary_t::DISPLACEMENT || bounds_cond[b].type[1] == boundary_t::DISPLACEMENT) {
                    bool push = true;
                    for(const std::vector<I>& bound : kinematic_nodes)
                        push = push && std::find(bound.cbegin(), bound.cend(), mesh().node_number(b, el, i)) == bound.cend();
                    if(push)
                        kinematic_nodes[b].push_back(mesh().node_number(b, el, i));
                }
            });

        for(size_t b = 0; b < kinematic_nodes.size(); ++b)
            for(const I node : kinematic_nodes[b]) 
                for(size_t comp = 0; comp < 2; ++comp)
                    if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT) {
                        const T temp = bounds_cond[b].func[comp](mesh().node(node));
                        for(typename Eigen::SparseMatrix<T>::InnerIterator it(K_bound, 2*node+comp); it; ++it)
                            f[it.row()] -= temp * it.value();
                    }

        // Повторный проход для корректировки
        for(size_t b = 0; b < kinematic_nodes.size(); ++b)
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT)
                    for(const I node : kinematic_nodes[b])
                        f[2*node+comp] = bounds_cond[b].func[comp](mesh().node(node));
    }

    template<class Vector>
    std::array<std::vector<std::array<T, 3>>, 2> strains_and_stress_loc(const Vector& displacement) const {
        std::vector<uint8_t> repeating(mesh().nodes_count(), 0); // Подсчёт повторений узла.
                                                                 // Будем надеяться, что один узел может принадлежать не более 255 элементам
        std::vector<std::array<T, 3>> strain(mesh().nodes_count(), std::array<T, 3>{}),
                                      stress(mesh().nodes_count(), std::array<T, 3>{});
        for(size_t el = 0; el < mesh().elements_count(); ++el) {
            const auto& e = mesh().element_2d(el);
            for(size_t i = 0; i < e->nodes_count(); ++i) {
                ++repeating[mesh().node_number(el, i)];
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
            strain[i][_11] /=   repeating[i];
            strain[i][_22] /=   repeating[i];
            strain[i][_12] /= 2*repeating[i];
            stress[i][_11] /=   repeating[i];
            stress[i][_22] /=   repeating[i];
            stress[i][_12] /= 2*repeating[i];
        }

        return {std::move(strain), std::move(stress)};
    }

    std::vector<std::array<T, 3>> approx_strain_in_quad(const std::vector<std::array<T, 3>>& strain) const {
        std::vector<std::array<T, 3>> strain_in_quad(quad_shift(mesh().elements_count()), std::array<T, 3>{});
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
        const T p2 = 1. - p1;
        const std::vector<std::array<T, 3>> strains_in_quad = approx_strain_in_quad(strains);
        for(size_t node = 0; node < mesh().nodes_count(); ++node)
            for(const auto elNL : mesh().node_neighbors(node)) {
                const auto& eNL = mesh().element_2d(mesh().element_2d_type(elNL));
                for(size_t q = 0, shift = quad_shift(elNL); q < eNL->qnodes_count(); ++q, ++shift) {
                    const T influence_weight = p2 * eNL->weight(q) * jacobian(shift) * influence_fun(quad_coord(shift), mesh().node(node));
                    stress[node][0] += influence_weight * (_D[0] * strains_in_quad[shift][0] + _D[1] * strains_in_quad[shift][1]);
                    stress[node][1] += influence_weight * (_D[1] * strains_in_quad[shift][0] + _D[0] * strains_in_quad[shift][1]);
                    stress[node][2] += influence_weight *  _D[2] * strains_in_quad[shift][2];
                }
            }
    }

public:
    explicit structural_solver(const mesh::mesh_2d<T, I>& mesh, const parameters<T>& params) :
        _base{mesh},
        _D{hooke_matrix(params)} {}

    explicit structural_solver(mesh::mesh_2d<T, I>&& mesh, const parameters<T>& params) :
        _base{std::move(mesh)},
        _D{hooke_matrix(params)} {}

    ~structural_solver() override = default;

    template<class Vector>
    void save_as_vtk(const std::string& path, const Vector& U,
                     const std::vector<std::array<T, 3>>& strain, const std::vector<std::array<T, 3>>& stress) const;

    template<class Influence_Function>
    Eigen::Matrix<T, Eigen::Dynamic, 1> stationary(
        const std::vector<boundary_condition<T>> &bounds_cond, //const distributed_load<T>& right_part,
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
        fout << sqrt(stress[i][_11] * stress[i][_11] + 
                     stress[i][_22] * stress[i][_22] -
                     stress[i][_11] * stress[i][_22] + 
                 3 * stress[i][_12] * stress[i][_12]) << '\n';
}

template<class T, class I>
template<class Influence_Function>
Eigen::Matrix<T, Eigen::Dynamic, 1> structural_solver<T, I>::stationary(
    const std::vector<boundary_condition<T>> &bounds_cond, //const distributed_load<T>& right_part,
    const T p1, const Influence_Function& influence_fun) {
    double time = omp_get_wtime();
    Eigen::SparseMatrix<T, Eigen::ColMajor, I> K      (2*mesh().nodes_count(), 2*mesh().nodes_count()),
                                               K_bound(2*mesh().nodes_count(), 2*mesh().nodes_count());
    create_matrix(K, K_bound, bounds_cond, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    //integrate_right_part(mesh, right_part, f);

    time = omp_get_wtime();
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(2*mesh().nodes_count());
    integrate_boundary_pressure(f, bounds_cond);
    displacement_on_boundary(f, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::PardisoLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;
    
    return std::move(u);
}

template<class T, class I>
template<class Vector, class Influence_Function>
std::array<std::vector<std::array<T, 3>>, 2> structural_solver<T, I>::strains_and_stress(
    const Vector& displacement, const T p1, const Influence_Function& influence_fun) const {
    auto [strain, stress] = strains_and_stress_loc(displacement);
    if(p1 < MAX_LOCAL_WEIGHT) { // Нелокальная задача
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
        //std::cout << el << " " << std::endl;
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