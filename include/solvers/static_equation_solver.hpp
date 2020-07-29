#ifndef STATIC_EQUATION_SOLVER_HPP
#define STATIC_EQUATION_SOLVER_HPP

#include <tuple>
#include <functional>
#include <algorithm>
#include "omp.h"
#include "finite_element_routine.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/PardisoSupport"

namespace nonlocal::structural {

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
        func = { [](const std::array<T, 2>&) noexcept { return 0.; },
                 [](const std::array<T, 2>&) noexcept { return 0.; } };
    std::array<boundary_t, 2> type = { boundary_t::PRESSURE, boundary_t::PRESSURE };
};

template<class T>
struct distributed_load {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    std::array<std::function<T(const std::array<T, 2>&)>, 2> 
        func = { [](const std::array<T, 2>&) noexcept { return 0.; },
                 [](const std::array<T, 2>&) noexcept { return 0.; } };
};

class _structural : protected _finite_element_routine {
protected:
    explicit _structural() noexcept = default;

    // Матрица Гука, которая имеет следующий портрет:
    // arr[0] arr[1]   0
    // arr[1] arr[0]   0
    //   0      0    arr[2]
    template<class T>
    static std::array<T, 3> hooke_matrix(const parameters<T>& params) noexcept {
        return {             params.E / (T(1) - params.nu*params.nu),
                 params.nu * params.E / (T(1) - params.nu*params.nu),
                    T(0.5) * params.E / (T(1) + params.nu) };
    }

    template<bool Proj, bool Form, class T, class Finite_Element_2D_Ptr>
    static T integrate_loc(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, 
                           const std::vector<std::array<T, 4>>& jacobi_matrices, size_t quad_shift,
                           const std::array<T, 3>& D) {
        T integral = 0;
        static constexpr size_t k = Proj ^ Form;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += (D[k] * dNd< Proj>(e, i, q, jacobi_matrices[quad_shift]) * dNd< Form>(e, j, q, jacobi_matrices[quad_shift]) +
                         D[2] * dNd<!Proj>(e, i, q, jacobi_matrices[quad_shift]) * dNd<!Form>(e, j, q, jacobi_matrices[quad_shift])) *
                        e->weight(q) / jacobian(jacobi_matrices[quad_shift]);
        return integral;
    }

    template<bool Proj, bool Form, class T, class Finite_Element_2D_Ptr, class Influence_Function>
    static T integrate_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,
                              const Finite_Element_2D_Ptr& eNL, const size_t jNL,
                              const std::vector<std::array<T, 2>>& quad_coords,
                              const std::vector<std::array<T, 4>>& jacobi_matrices,
                              size_t shiftL, size_t shiftNL,
                              const Influence_Function& influence_function,
                              const std::array<T, 3> &D) {
        T integral = 0;
        const size_t sub_shift = shiftNL;
        static constexpr size_t k = Proj ^ Form;
        for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL) {
            T int_with_weight_x = 0, int_with_weight_y = 0;
            for(size_t qNL = 0, shiftNL = sub_shift; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL) {
                const T influence_weight = eNL->weight(qNL) * influence_function(quad_coords[shiftL], quad_coords[shiftNL]);
                int_with_weight_x += influence_weight * dNd< Form>(eNL, jNL, qNL, jacobi_matrices[shiftNL]);
                int_with_weight_y += influence_weight * dNd<!Form>(eNL, jNL, qNL, jacobi_matrices[shiftNL]);
            }
            integral += eL->weight(qL) *
                        (D[k] * int_with_weight_x * dNd< Proj>(eL, iL, qL, jacobi_matrices[shiftL]) +
                         D[2] * int_with_weight_y * dNd<!Proj>(eL, iL, qL, jacobi_matrices[shiftL]));
        }
        return integral;
    }

    template<class Type, class Index>
    static std::array<std::vector<Index>, 4>
        mesh_analysis(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<bool>& inner_nodes, const bool nonlocal) {
        std::vector<Index> shifts_loc      (mesh.elements_count()+1, 0), 
                           shifts_bound_loc(mesh.elements_count()+1, 0),
                           shifts_nonloc, shifts_bound_nonloc;

        const auto counter_loc = 
            [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc]
            (size_t i, size_t j, size_t el, component proj, component form) {
                size_t row = 2 * mesh.node_number(el, i) + size_t(proj),
                       col = 2 * mesh.node_number(el, j) + size_t(form);
                if(row >= col) {
                    if(inner_nodes[row] && inner_nodes[col])
                        ++shifts_loc[el+1];
                    else if(row != col)
                        ++shifts_bound_loc[el+1];
                }
            };

        mesh_run_loc(mesh, 
            [&counter_loc](size_t i, size_t j, size_t el) {
                counter_loc(i, j, el, component::X, component::X);
                counter_loc(i, j, el, component::X, component::Y);
                counter_loc(i, j, el, component::Y, component::X);
                counter_loc(i, j, el, component::Y, component::Y);
            });

        shifts_loc[0] = std::count(inner_nodes.cbegin(), inner_nodes.cend(), false);
        for(size_t i = 1; i < shifts_loc.size(); ++i) {
            shifts_loc[i] += shifts_loc[i-1];
            shifts_bound_loc[i] += shifts_bound_loc[i-1];
        }

        if(nonlocal) {
            shifts_nonloc.resize(mesh.elements_count()+1, 0);
            shifts_bound_nonloc.resize(mesh.elements_count()+1, 0);

            const auto counter_nonloc =
                [&mesh, &inner_nodes, &shifts_nonloc, &shifts_bound_nonloc]
                (size_t iL, size_t jNL, size_t elL, size_t elNL, component proj, component form) {
                    size_t row = 2 * mesh.node_number(elL , iL ) + size_t(proj),
                           col = 2 * mesh.node_number(elNL, jNL) + size_t(form);
                    if(row >= col) {
                        if(inner_nodes[row] && inner_nodes[col])
                            ++shifts_nonloc[elL+1];
                        else if(row != col)
                            ++shifts_bound_nonloc[elL+1];
                    }
                };

            mesh_run_nonloc(mesh,
                [&counter_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL) {
                    counter_nonloc(iL, jNL, elL, elNL, component::X, component::X);
                    counter_nonloc(iL, jNL, elL, elNL, component::X, component::Y);
                    counter_nonloc(iL, jNL, elL, elNL, component::Y, component::X);
                    counter_nonloc(iL, jNL, elL, elNL, component::Y, component::Y);
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

    template<class Type, class Index, class Influence_Function>
    static std::array<std::vector<Eigen::Triplet<Type, Index>>, 2>
        triplets_fill(const mesh::mesh_2d<Type, Index>& mesh,
                      const std::vector<Index>& shifts_quad,
                      const std::vector<std::array<Type, 2>>& all_quad_coords,
                      const std::vector<std::array<Type, 4>>& all_jacobi_matrices,
                      const std::vector<boundary_condition<Type>> &bounds_cond,
                      const parameters<Type>& params,
                      const Type p1, const Influence_Function& influence_fun) {
        const bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
        std::vector<bool> inner_nodes(2*mesh.nodes_count(), true);
        boundary_nodes_run(mesh, [&mesh, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT)
                    inner_nodes[2*mesh.node_number(b, el, i)+comp] = false;
        });

        auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, inner_nodes, nonlocal);
        std::vector<Eigen::Triplet<Type, Index>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                                 triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
        for(size_t i = 0, j = 0; i < inner_nodes.size(); ++i)
            if(!inner_nodes[i])
                triplets[j++] = Eigen::Triplet<Type, Index>(i, i, 1.);

        const std::array<Type, 3> D = hooke_matrix(params);
        const auto filler_loc =
            [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc, &triplets, &triplets_bound, &shifts_quad, &all_jacobi_matrices, p1, &D]
            (size_t i, size_t j, size_t el, component proj, component form, const auto& integrate_rule) {
                size_t row = 2 * mesh.node_number(el, i) + size_t(proj),
                       col = 2 * mesh.node_number(el, j) + size_t(form);
                if(row >= col) {
                    Type integral = p1 * integrate_rule(mesh.element_2d(mesh.element_2d_type(el)), i, j, all_jacobi_matrices, shifts_quad[el], D);
                    if(inner_nodes[row] && inner_nodes[col])
                        triplets[shifts_loc[el]++] = Eigen::Triplet<Type, Index>(row, col, integral);
                    else if(row != col)
                        triplets_bound[shifts_bound_loc[el]++] = inner_nodes[col] ? Eigen::Triplet<Type, Index>(col, row, integral) :
                                                                                    Eigen::Triplet<Type, Index>(row, col, integral);
                }
            };

        mesh_run_loc(mesh, 
            [&filler_loc](size_t i, size_t j, size_t el) {
                using Finite_Element_2D_Ptr = typename mesh::mesh_2d<Type, Index>::fe_2d_ptr;
                filler_loc(i, j, el, X, X, integrate_loc<X, X, Type, Finite_Element_2D_Ptr>);
                filler_loc(i, j, el, X, Y, integrate_loc<X, Y, Type, Finite_Element_2D_Ptr>);
                filler_loc(i, j, el, Y, X, integrate_loc<Y, X, Type, Finite_Element_2D_Ptr>);
                filler_loc(i, j, el, Y, Y, integrate_loc<Y, Y, Type, Finite_Element_2D_Ptr>);
            });

        if(nonlocal) {
            const auto filler_nonloc =
                [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
                &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1, &D]
                (size_t iL, size_t jNL, size_t elL, size_t elNL, component proj, component form, const auto& integrate_rule) {
                    size_t row = 2 * mesh.node_number(elL,  iL ) + size_t(proj),
                           col = 2 * mesh.node_number(elNL, jNL) + size_t(form);
                    if(row >= col) {
                        Type integral = p2 * integrate_rule(mesh.element_2d(mesh.element_2d_type(elL )), iL,
                                                            mesh.element_2d(mesh.element_2d_type(elNL)), jNL,
                                                            all_quad_coords, all_jacobi_matrices,
                                                            shifts_quad[elL], shifts_quad[elNL],
                                                            influence_fun, D);
                        if(inner_nodes[row] && inner_nodes[col])
                            triplets[shifts_nonloc[elL]++] = Eigen::Triplet<Type, Index>(row, col, integral);
                        else if(row != col)
                            triplets_bound[shifts_bound_nonloc[elL]++] = inner_nodes[col] ? Eigen::Triplet<Type, Index>(col, row, integral) :
                                                                                            Eigen::Triplet<Type, Index>(row, col, integral);
                    }
                };

            mesh_run_nonloc(mesh, 
                [&filler_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL) {
                    using Finite_Element_2D_Ptr = typename mesh::mesh_2d<Type, Index>::fe_2d_ptr;
                    filler_nonloc(iL, jNL, elL, elNL, X, X, integrate_nonloc<X, X, Type, Finite_Element_2D_Ptr, Influence_Function>);
                    filler_nonloc(iL, jNL, elL, elNL, X, Y, integrate_nonloc<X, Y, Type, Finite_Element_2D_Ptr, Influence_Function>);
                    filler_nonloc(iL, jNL, elL, elNL, Y, X, integrate_nonloc<Y, X, Type, Finite_Element_2D_Ptr, Influence_Function>);
                    filler_nonloc(iL, jNL, elL, elNL, Y, Y, integrate_nonloc<Y, Y, Type, Finite_Element_2D_Ptr, Influence_Function>);
                });
        }

        return {std::move(triplets), std::move(triplets_bound)};
    }

    template<class Type, class Index, class Influence_Function>
    static void create_matrix(Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K, Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K_bound,
                              const mesh::mesh_2d<Type, Index>& mesh, 
                              const std::vector<Index>& shifts_quad,
                              const std::vector<std::array<Type, 2>>& all_quad_coords,
                              const std::vector<std::array<Type, 4>>& all_jacobi_matrices,
                              const std::vector<boundary_condition<Type>>& bounds_cond, const parameters<Type>& params,
                              const Type p1, const Influence_Function& influence_fun) {
        double time = omp_get_wtime();
        auto [triplets, triplets_bound] = triplets_fill(mesh, shifts_quad, all_quad_coords, all_jacobi_matrices,
                                                        bounds_cond, params, p1, influence_fun);
        std::cout << omp_get_wtime() - time << std::endl
                  << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;

        K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
        triplets_bound.reserve(0);
        K.setFromTriplets(triplets.cbegin(), triplets.cend());
        std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
    }

    template<class Type, class Index, class Vector>
    static void integrate_boundary_pressure(Vector& f, const mesh::mesh_2d<Type, Index>& mesh,
                                            const std::vector<boundary_condition<Type>>& bounds_cond) {
        std::vector<std::array<Type, 2>> quad_nodes, jacobi_matrices;
        for(size_t b = 0; b < bounds_cond.size(); ++b)
            if(bounds_cond[b].type[0] == boundary_t::PRESSURE || bounds_cond[b].type[1] == boundary_t::PRESSURE)
                for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                    approx_quad_nodes_on_bound(quad_nodes, mesh, b, el);
                    approx_jacobi_matrices_on_bound(jacobi_matrices, mesh, b, el);
                    const auto& be = mesh.element_1d(mesh.element_1d_type(b, el));
                    for(size_t i = 0; i < be->nodes_count(); ++i) 
                        for(size_t comp = 0; comp < 2; ++comp)
                            if(bounds_cond[b].type[comp] == boundary_t::PRESSURE)
                                f[2*mesh.node_number(b, el, i)+comp] += 
                                    integrate_boundary_gradient(be, i, quad_nodes, jacobi_matrices, bounds_cond[b].func[comp]);
                }
    }

    // Учёт граничных условий первого рода.
    template<class Type, class Index, class Vector>
    static void displacement_on_boundary(Vector& f, const mesh::mesh_2d<Type, Index>& mesh,
                                         const std::vector<boundary_condition<Type>>& bounds_cond,
                                         const Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K_bound) {
        std::vector<std::vector<Index>> kinematic_nodes(mesh.boundary_groups_count());
        boundary_nodes_run(mesh, [&mesh, &bounds_cond, &kinematic_nodes](const size_t b, const size_t el, const size_t i) {
            if(bounds_cond[b].type[0] == boundary_t::DISPLACEMENT || bounds_cond[b].type[1] == boundary_t::DISPLACEMENT) {
                bool push = true;
                for(const std::vector<Index>& bound : kinematic_nodes)
                    push = push && std::find(bound.cbegin(), bound.cend(), mesh.node_number(b, el, i)) == bound.cend();
                if(push)
                    kinematic_nodes[b].push_back(mesh.node_number(b, el, i));
            }
        });

        for(size_t b = 0; b < kinematic_nodes.size(); ++b)
            for(const Index node : kinematic_nodes[b]) 
                for(size_t comp = 0; comp < 2; ++comp)
                    if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT) {
                        const Type temp = bounds_cond[b].func[comp](mesh.node(node));
                        for(typename Eigen::SparseMatrix<Type>::InnerIterator it(K_bound, 2*node+comp); it; ++it)
                            f[it.row()] -= temp * it.value();
                    }

        // Повторный проход для корректировки
        for(size_t b = 0; b < kinematic_nodes.size(); ++b)
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type[comp] == boundary_t::DISPLACEMENT)
                    for(const Index node : kinematic_nodes[b])
                        f[2*node+comp] = bounds_cond[b].func[comp](mesh.node(node));
    }

    template<class Type, class Index, class Vector>
    static std::array<std::vector<std::array<Type, 3>>, 2>
        strains_and_stress_loc(const mesh::mesh_2d<Type, Index>& mesh, const std::array<Type, 3>& D, const Vector& displacement) {
        std::vector<uint8_t> repeating(mesh.nodes_count(), 0); // Подсчёт повторений узла.
                                                               // Будем надеяться, что один узел может принадлежать не более 255 элементам
        std::vector<std::array<Type, 3>> strain(mesh.nodes_count(), std::array<Type, 3>{}),
                                         stress(mesh.nodes_count(), std::array<Type, 3>{});
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t i = 0; i < e->nodes_count(); ++i) {
                ++repeating[mesh.node_number(el, i)];
                std::array<Type, 4> jacobi_matrix = {};
                for(size_t j = 0; j < e->nodes_count(); ++j) {
                    jacobi_matrix[0] += mesh.node(mesh.node_number(el, j))[0] * e->Nxi (j, e->node(i));
                    jacobi_matrix[1] += mesh.node(mesh.node_number(el, j))[0] * e->Neta(j, e->node(i));
                    jacobi_matrix[2] += mesh.node(mesh.node_number(el, j))[1] * e->Nxi (j, e->node(i));
                    jacobi_matrix[3] += mesh.node(mesh.node_number(el, j))[1] * e->Neta(j, e->node(i));
                }

                std::array<Type, 3> strain_loc = {};
                for(size_t j = 0; j < e->nodes_count(); ++j) {
                    const Type jac = jacobian(jacobi_matrix),
                               dx1 =  jacobi_matrix[3] * e->Nxi(j, e->node(i)) - jacobi_matrix[2] * e->Neta(j, e->node(i)),
                               dx2 = -jacobi_matrix[1] * e->Nxi(j, e->node(i)) + jacobi_matrix[0] * e->Neta(j, e->node(i));
                    strain_loc[0] +=  dx1 * displacement[2*mesh.node_number(el, j)  ]  / jac;
                    strain_loc[1] +=  dx2 * displacement[2*mesh.node_number(el, j)+1]  / jac;
                    strain_loc[2] += (dx2 * displacement[2*mesh.node_number(el, j)  ] +
                                      dx1 * displacement[2*mesh.node_number(el, j)+1]) / jac;
                }

                stress[mesh.node_number(el, i)][0] += D[0] * strain_loc[0] + D[1] * strain_loc[1];
                stress[mesh.node_number(el, i)][1] += D[1] * strain_loc[0] + D[0] * strain_loc[1];
                stress[mesh.node_number(el, i)][2] += D[2] * strain_loc[2];
                for(size_t j = 0; j < 3; ++j)
                    strain[mesh.node_number(el, i)][j] += strain_loc[j];
            }
        }

        for(size_t i = 0; i < mesh.nodes_count(); ++i) {
            strain[i][0] /=   repeating[i];
            strain[i][1] /=   repeating[i];
            strain[i][2] /= 2*repeating[i];
            stress[i][0] /=   repeating[i];
            stress[i][1] /=   repeating[i];
            stress[i][2] /= 2*repeating[i];
        }

        return {std::move(strain), std::move(stress)};
    }

    template<class Type, class Index>
    static std::vector<std::array<Type, 3>>
        approx_strains_in_quad(const mesh::mesh_2d<Type, Index> &mesh, const std::vector<Index> &shifts,
                               const std::vector<std::array<Type, 3>>& strains) {
        std::vector<std::array<Type, 3>> strains_in_quad(shifts.back(), std::array<Type, 3>{});
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t q = 0, shift = shifts[el]; q < e->qnodes_count(); ++q, ++shift)
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 3; ++comp)
                        strains_in_quad[shift][comp] += strains[mesh.node_number(el, i)][comp] * e->qN(i, q);
        }
        return std::move(strains_in_quad);
    }

    template<class Type, class Index, class Influence_Function>
    static void stress_nonloc(std::vector<std::array<Type, 3>>& stress,
                              const std::vector<std::array<Type, 3>>& strain,
                              const mesh::mesh_2d<Type, Index>& mesh, const std::array<Type, 3>& D,
                              const Type p1, const Influence_Function& influence_fun) {
        const Type p2 = 1. - p1;
        const std::vector<Index> shifts_quad = quadrature_shifts_init(mesh);
        const std::vector<std::array<Type, 2>> all_quad_coords = approx_all_quad_nodes(mesh, shifts_quad);
        const std::vector<std::array<Type, 3>> strains_in_quad = approx_strains_in_quad(mesh, shifts_quad, strain);
        const std::vector<std::array<Type, 4>> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);
        for(size_t node = 0; node < mesh.nodes_count(); ++node)
            for(const auto elNL : mesh.node_neighbors(node)) {
                const auto& eNL = mesh.element_2d(mesh.element_2d_type(elNL));
                for(size_t q = 0, shift = shifts_quad[elNL]; q < eNL->qnodes_count(); ++q, ++shift) {
                    const Type influence_weight = p2 * eNL->weight(q) * jacobian(all_jacobi_matrices[shift]) *
                                                  influence_fun(all_quad_coords[shift], mesh.node(node));
                    stress[node][0] += influence_weight * (D[0] * strains_in_quad[shift][0] + D[1] * strains_in_quad[shift][1]);
                    stress[node][1] += influence_weight * (D[1] * strains_in_quad[shift][0] + D[0] * strains_in_quad[shift][1]);
                    stress[node][2] += influence_weight *  D[2] * strains_in_quad[shift][2];
                }
            }
    }

public:
    template<class Type, class Index, class Vector>
    friend void save_as_vtk(const std::string& path, const mesh::mesh_2d<Type, Index>& mesh, const Vector& U,
                            const std::vector<std::array<Type, 3>>& strain,
                            const std::vector<std::array<Type, 3>>& stress);

    template<class Type, class Index, class Influence_Function>
    friend Eigen::Matrix<Type, Eigen::Dynamic, 1>
        stationary(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>> &bounds_cond,
                   const parameters<Type>& params, //const distributed_load<Type>& right_part,
                   const Type p1, const Influence_Function& influence_fun);

    template<class Type, class Index, class Vector, class Influence_Function>
    friend std::array<std::vector<std::array<Type, 3>>, 2>
        strains_and_stress(const mesh::mesh_2d<Type, Index>& mesh, const parameters<Type>& params, const Vector& displacement,
                           const Type p1, const Influence_Function& influence_fun);
};

template<class Type, class Index, class Vector>
void save_as_vtk(const std::string& path, const mesh::mesh_2d<Type, Index>& mesh, const Vector& U,
                 const std::vector<std::array<Type, 3>>& strain,
                 const std::vector<std::array<Type, 3>>& stress) {
    static constexpr std::string_view data_type = std::is_same_v<Type, double> ? "double" : "float";

    if(2 * mesh.nodes_count() != size_t(U.size()))
        throw std::domain_error{"2 * mesh.nodes_count() != U.size()."};

    std::ofstream fout{path};
    fout.precision(20);

    mesh.save_as_vtk(fout);

    fout << "POINT_DATA " << mesh.nodes_count() << std::endl;
    fout << "VECTORS Displacement " << data_type << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << U[2*i] << ' ' << U[2*i+1] << " 0\n";

    static constexpr std::array<std::string_view, 3>
        strain_number = {"strain11", "strain22", "strain12"},
        stress_number = {"stress11", "stress22", "stress12"};

    for(size_t comp = 0; comp < 3; ++comp) {
        fout << "SCALARS " << strain_number[comp] << ' ' << data_type << " 1" << std::endl
             << "LOOKUP_TABLE default" << std::endl;
        for(size_t i = 0; i < mesh.nodes_count(); ++i)
            fout << strain[i][comp] << '\n';
    }

    for(size_t comp = 0; comp < 3; ++comp) {
        fout << "SCALARS " << stress_number[comp] << ' ' << data_type << " 1" << std::endl
             << "LOOKUP_TABLE default" << std::endl;
        for(size_t i = 0; i < mesh.nodes_count(); ++i)
            fout << stress[i][comp] << '\n';
    }
}

template<class Type, class Index, class Influence_Function>
Eigen::Matrix<Type, Eigen::Dynamic, 1>
    stationary(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>> &bounds_cond,
               const parameters<Type>& params, //const distributed_load<Type>& right_part,
               const Type p1, const Influence_Function& influence_fun) {
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                                      K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(2*mesh.nodes_count());
    
    std::cout << "Quadratures data init: ";
    double time = omp_get_wtime();
    const std::vector<Index> shifts_quad = _structural::quadrature_shifts_init(mesh);
    const std::vector<std::array<Type, 2>> all_quad_coords = _structural::approx_all_quad_nodes(mesh, shifts_quad);
    const std::vector<std::array<Type, 4>> all_jacobi_matrices = _structural::approx_all_jacobi_matrices(mesh, shifts_quad);
    std::cout << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    _structural::create_matrix(K, K_bound, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices,
                               bounds_cond, params, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    //integrate_right_part(mesh, right_part, f);

    time = omp_get_wtime();
    _structural::integrate_boundary_pressure(f, mesh, bounds_cond);
    _structural::displacement_on_boundary(f, mesh, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;
    
    return u;
}

template<class Type, class Index, class Vector, class Influence_Function>
std::array<std::vector<std::array<Type, 3>>, 2>
    strains_and_stress(const mesh::mesh_2d<Type, Index>& mesh, const parameters<Type>& params, const Vector& displacement,
                       const Type p1, const Influence_Function& influence_fun) {
    const std::array<Type, 3> D = _structural::hooke_matrix(params);
    auto [strain, stress] = _structural::strains_and_stress_loc(mesh, D, displacement);

    if(p1 < _structural::MAX_LOCAL_WEIGHT) { // Нелокальная задача
        for(size_t i = 0; i < mesh.nodes_count(); ++i)
            for(size_t j = 0; j < 3; ++j)
                stress[i][j] *= p1;
        _structural::stress_nonloc(stress, strain, mesh, D, p1, influence_fun);
    }

    return {std::move(strain), std::move(stress)};
}

}

#endif