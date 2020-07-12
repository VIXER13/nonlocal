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

namespace static_equation_with_nonloc
{

enum class boundary_type : uint8_t {DISPLACEMENT, PRESSURE};

template<class Type>
struct boundary_condition
{
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    std::function<Type(Type, Type)> func_x = [](Type, Type) { return 0.; },
                                    func_y = [](Type, Type) { return 0.; };
    boundary_type type_x = boundary_type::PRESSURE,
                  type_y = boundary_type::PRESSURE;
};

template<class Type>
struct distributed_load
{
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    std::function<Type(Type, Type)> func_x = [](Type, Type) { return 0.; },
                                    func_y = [](Type, Type) { return 0.; };
};

template<class Type>
struct parameters
{
    Type nu = 0, // Коэффициент Пуассона
         E  = 0;  // Модуль Юнга
};

enum class component : uint8_t {X, Y};

template<class Type>
static std::array<Type, 3> hooke_matrix(const parameters<Type> &params)
{
    return {            params.E / (Type(1.) - params.nu*params.nu),
            params.nu * params.E / (Type(1.) - params.nu*params.nu),
                      Type(0.5) * params.E / (Type(1.) + params.nu) };
}

template<class Type, component Projection, component Form>
static Type integrate_loc(const finite_element::element_2d_integrate_base<Type> *const e,
                          const size_t i, const size_t j, const matrix<Type> &jacobi_matrices, size_t shift,
                          const std::array<Type, 3> &D)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
    {
        if constexpr (Projection == component::X && Form == component::X)
            integral += (D[0] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2)) +
                         D[2] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0))) /
                        (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);

        if constexpr (Projection == component::X && Form == component::Y)
            integral += (D[1] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0)) +
                         D[2] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2))) /
                        (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);

        if constexpr (Projection == component::Y && Form == component::X)
            integral += (D[1] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2)) +
                         D[2] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0))) /
                        (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);

        if constexpr (Projection == component::Y && Form == component::Y)
            integral += (D[0] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0)) +
                         D[2] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2))) /
                        (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
    }
    return integral;
}

template<class Type, component Projection, component Form>
static Type integrate_nonloc(const finite_element::element_2d_integrate_base<Type> *const eL,
                             const finite_element::element_2d_integrate_base<Type> *const eNL,
                             const size_t iL, const size_t jNL, size_t shiftL, size_t shiftNL,
                             const matrix<Type> &coords, const matrix<Type> &jacobi_matrices,
                             const std::function<Type(Type, Type, Type, Type)> &influence_fun,
                             const std::array<Type, 3> &D)
{
    const size_t sub = shiftNL;
    Type integral = 0., int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
    for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL)
    {
        int_with_weight_x = int_with_weight_y = 0.;
        for(size_t qNL = 0, shiftNL = sub; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL)
        {
            finit = eNL->weight(qNL) * influence_fun(coords(shiftL, 0), coords(shiftNL, 0), coords(shiftL, 1), coords(shiftNL, 1));
            if constexpr (Form == component::X)
            {
                int_with_weight_x += finit * ( eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 3) - eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 2));
                int_with_weight_y += finit * (-eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 1) + eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 0));
            }
            else
            {
                int_with_weight_x += finit * (-eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 1) + eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 0));
                int_with_weight_y += finit * ( eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 3) - eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 2));
            }
        }

        if constexpr (Projection == component::X)
            integral += eL->weight(qL) *
                        (D[size_t(Form)] * int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)) +
                         D[2]            * int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)));
        else
            integral += eL->weight(qL) *
                        (D[1-size_t(Form)] * int_with_weight_x * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)) +
                         D[2]              * int_with_weight_y * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)));
    }
    return integral;
}

template<class Type, class Index>
static std::vector<bool> inner_nodes_vector(const mesh_2d<Type, Index> &mesh, const std::vector<boundary_condition<Type>> &bounds_cond)
{
    std::vector<bool> inner_nodes(2*mesh.nodes_count(), true);
    for(size_t b = 0; b < bounds_cond.size(); ++b)
    {
        if(bounds_cond[b].type_x == boundary_type::DISPLACEMENT)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                inner_nodes[2*(*node)] = false;
                
        if(bounds_cond[b].type_y == boundary_type::DISPLACEMENT)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                inner_nodes[2*(*node)+1] = false;
    }
    return std::move(inner_nodes);
}

template<class Type, class Index>
static std::vector<std::vector<Index>>
    kinematic_nodes_vectors(const mesh_2d<Type, Index> &mesh, const std::vector<boundary_condition<Type>> &bounds_cond)
{
    std::vector<std::vector<Index>> kinematic_nodes(bounds_cond.size());
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(bounds_cond[b].type_x == boundary_type::DISPLACEMENT || bounds_cond[b].type_y == boundary_type::DISPLACEMENT)
            for(auto [node, k] = std::make_tuple(mesh.boundary(b).cbegin(), size_t(0)); node != mesh.boundary(b).cend(); ++node)
            {
                for(k = 0; k < kinematic_nodes.size(); ++k)
                    if(std::find(kinematic_nodes[k].cbegin(), kinematic_nodes[k].cend(), *node) != kinematic_nodes[k].cend())
                        break;
                if(k == kinematic_nodes.size())
                    kinematic_nodes[b].push_back(*node);
            }
    return std::move(kinematic_nodes);
}

template<class Type, class Index>
static std::array<std::vector<Index>, 4>
    mesh_analysis(const mesh_2d<Type, Index> &mesh, const std::vector<bool> &inner_nodes, const bool nonlocal)
{
    std::vector<Index> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0),
                       shifts_nonloc, shifts_bound_nonloc;

    const auto counter_loc = 
        [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc]
        (size_t i, size_t j, size_t el, component projection, component form)
        {
            size_t row = 2 * mesh.node_number(el, i) + size_t(projection),
                   col = 2 * mesh.node_number(el, j) + size_t(form);
            if(row >= col)
            {
                if(inner_nodes[row] && inner_nodes[col])
                    ++shifts_loc[el+1];
                else if(row != col)
                    ++shifts_bound_loc[el+1];
            }
        };

    mesh_run_loc(mesh, [&counter_loc](size_t i, size_t j, size_t el)
                       {
                           counter_loc(i, j, el, component::X, component::X);
                           counter_loc(i, j, el, component::X, component::Y);
                           counter_loc(i, j, el, component::Y, component::X);
                           counter_loc(i, j, el, component::Y, component::Y);
                       });

    shifts_loc[0] = std::count(inner_nodes.cbegin(), inner_nodes.cend(), false);
    for(size_t i = 1; i < shifts_loc.size(); ++i)
    {
        shifts_loc[i] += shifts_loc[i-1];
        shifts_bound_loc[i] += shifts_bound_loc[i-1];
    }

    if(nonlocal)
    {
        shifts_nonloc.resize(mesh.elements_count()+1, 0);
        shifts_bound_nonloc.resize(mesh.elements_count()+1, 0);

        const auto counter_nonloc =
            [&mesh, &inner_nodes, &shifts_nonloc, &shifts_bound_nonloc]
            (size_t iL, size_t jNL, size_t elL, size_t elNL, component projection, component form)
            {
                size_t row = 2 * mesh.node_number(elL , iL ) + size_t(projection),
                       col = 2 * mesh.node_number(elNL, jNL) + size_t(form);
                if(row >= col)
                {
                    if(inner_nodes[row] && inner_nodes[col])
                        ++shifts_nonloc[elL+1];
                    else if(row != col)
                        ++shifts_bound_nonloc[elL+1];
                }
            };

        mesh_run_nonloc(mesh, [&counter_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
                              {
                                  counter_nonloc(iL, jNL, elL, elNL, component::X, component::X);
                                  counter_nonloc(iL, jNL, elL, elNL, component::X, component::Y);
                                  counter_nonloc(iL, jNL, elL, elNL, component::Y, component::X);
                                  counter_nonloc(iL, jNL, elL, elNL, component::Y, component::Y);
                              });

        shifts_nonloc[0] = shifts_loc.back();
        shifts_bound_nonloc[0] = shifts_bound_loc.back();
        for(size_t i = 1; i < shifts_nonloc.size(); ++i)
        {
            shifts_nonloc[i] += shifts_nonloc[i-1];
            shifts_bound_nonloc[i] += shifts_bound_nonloc[i-1];
        }
    }

    return {std::move(shifts_loc), std::move(shifts_bound_loc), std::move(shifts_nonloc), std::move(shifts_bound_nonloc)};
}

template<class Type, class Index>
static std::array<std::vector<Eigen::Triplet<Type, Index>>, 2>
    triplets_fill(const mesh_2d<Type, Index> &mesh, const std::vector<bool> &inner_nodes, const parameters<Type> &params,
                  const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    static constexpr Type MAX_LOCAL_WEIGHT = 0.999;
    bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
    auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, inner_nodes, nonlocal);
    std::vector<Eigen::Triplet<Type, Index>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                             triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
    std::cout << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;
    for(size_t i = 0, j = 0; i < inner_nodes.size(); ++i)
        if(!inner_nodes[i])
            triplets[j++] = Eigen::Triplet<Type, Index>(i, i, 1.);

    const std::vector<Index> shifts_quad = quadrature_shifts_init(mesh);
    const matrix<Type> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);
    const std::array<Type, 3> D = hooke_matrix(params);

    const auto filler_loc =
        [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc, &triplets, &triplets_bound, &shifts_quad, &all_jacobi_matrices, p1, &D]
        (size_t i, size_t j, size_t el, component projection, component form, const auto& integrate_rule)
        {
            size_t row = 2 * mesh.node_number(el, i) + size_t(projection),
                   col = 2 * mesh.node_number(el, j) + size_t(form);
            if(row >= col)
            {
                Type integral = p1 * integrate_rule(mesh.element_2d(mesh.element_type(el)), i, j, all_jacobi_matrices, shifts_quad[el], D);
                if(inner_nodes[row] && inner_nodes[col])
                    triplets[shifts_loc[el]++] = Eigen::Triplet<Type, Index>(row, col, integral);
                else if(row != col)
                    triplets_bound[shifts_bound_loc[el]++] = inner_nodes[col] ? Eigen::Triplet<Type, Index>(col, row, integral) :
                                                                                Eigen::Triplet<Type, Index>(row, col, integral);
            }
        };

    mesh_run_loc(mesh, [&filler_loc](size_t i, size_t j, size_t el)
                       {
                            filler_loc(i, j, el, component::X, component::X, integrate_loc<Type, component::X, component::X>);
                            filler_loc(i, j, el, component::X, component::Y, integrate_loc<Type, component::X, component::Y>);
                            filler_loc(i, j, el, component::Y, component::X, integrate_loc<Type, component::Y, component::X>);
                            filler_loc(i, j, el, component::Y, component::Y, integrate_loc<Type, component::Y, component::Y>);
                       });

    if(nonlocal)
    {
        const matrix<Type> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts_quad);

        const auto filler_nonloc =
            [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
             &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1, &D]
            (size_t iL, size_t jNL, size_t elL, size_t elNL, component projection, component form, const auto& integrate_rule)
            {
                size_t row = 2 * mesh.node_number(elL,  iL ) + size_t(projection),
                       col = 2 * mesh.node_number(elNL, jNL) + size_t(form);
                if(row >= col)
                {
                    Type integral = p2 * integrate_rule(mesh.element_2d(mesh.element_type(elL )),
                                                        mesh.element_2d(mesh.element_type(elNL)), 
                                                        iL, jNL, shifts_quad[elL], shifts_quad[elNL],
                                                        all_quad_coords, all_jacobi_matrices, influence_fun, D);
                    if(inner_nodes[row] && inner_nodes[col])
                        triplets[shifts_nonloc[elL]++] = Eigen::Triplet<Type, Index>(row, col, integral);
                    else if(row != col)
                        triplets_bound[shifts_bound_nonloc[elL]++] = inner_nodes[col] ? Eigen::Triplet<Type, Index>(col, row, integral) :
                                                                                        Eigen::Triplet<Type, Index>(row, col, integral);
                }
            };

        mesh_run_nonloc(mesh, [&filler_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
                              {
                                  filler_nonloc(iL, jNL, elL, elNL, component::X, component::X, integrate_nonloc<Type, component::X, component::X>);
                                  filler_nonloc(iL, jNL, elL, elNL, component::X, component::Y, integrate_nonloc<Type, component::X, component::Y>);
                                  filler_nonloc(iL, jNL, elL, elNL, component::Y, component::X, integrate_nonloc<Type, component::Y, component::X>);
                                  filler_nonloc(iL, jNL, elL, elNL, component::Y, component::Y, integrate_nonloc<Type, component::Y, component::Y>);
                              });
    }

    return {std::move(triplets), std::move(triplets_bound)};
}

template<int MatrixMajor, class Type, class Index>
static void create_matrix(const mesh_2d<Type, Index> &mesh, const parameters<Type> &params,
                          const std::vector<boundary_condition<Type>> &bounds_cond,
                          Eigen::SparseMatrix<Type, MatrixMajor, Index> &K, Eigen::SparseMatrix<Type, MatrixMajor, Index> &K_bound,
                          const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    double time = omp_get_wtime();
    auto [triplets, triplets_bound] = triplets_fill(mesh, inner_nodes_vector(mesh, bounds_cond), params, p1, influence_fun);
    std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;

    K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
    triplets_bound.clear();
    K.setFromTriplets(triplets.cbegin(), triplets.cend());
    std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
}

template<component Shift, class Type, class Index>
static void pressure_cond(const mesh_2d<Type, Index> &mesh, const size_t b, const std::function<Type(Type, Type)> &func,
                          Eigen::Matrix<Type, Eigen::Dynamic, 1> &f)
{
    matrix<Type> coords, jacobi_matrices;
    for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
    {
        const auto& be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
        approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
        approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
        for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
            f[2*mesh.boundary(b)(el, i)+Index(Shift)] += integrate_boundary_gradient(be, i, coords, jacobi_matrices, func);
    }
}

template<component Shift, int MatrixMajor, class Type, class Index>
static void displacement_cond(const mesh_2d<Type, Index> &mesh, const size_t b, const std::vector<std::vector<Index>> &temperature_nodes,
                              const std::function<Type(Type, Type)> &func,
                              const Eigen::SparseMatrix<Type, MatrixMajor, Index> &K_bound, Eigen::Matrix<Type, Eigen::Dynamic, 1> &f)
{
    for(auto node : temperature_nodes[b])
    {
        const Type temp = func(mesh.coord(node, 0), mesh.coord(node, 1));
        for(typename Eigen::SparseMatrix<Type, MatrixMajor, Index>::InnerIterator it(K_bound, 2*node+Index(Shift)); it; ++it)
            f[it.row()] -= temp * it.value();
    }
}

template<int MatrixMajor, class Type, class Index>
static void boundary_condition_calc(const mesh_2d<Type, Index> &mesh, const std::vector<std::vector<Index>> &temperature_nodes,
                                    const std::vector<boundary_condition<Type>> &bounds_cond,
                                    const Eigen::SparseMatrix<Type, MatrixMajor, Index> &K_bound,
                                    Eigen::Matrix<Type, Eigen::Dynamic, 1> &f)
{
    for(size_t b = 0; b < bounds_cond.size(); ++b)
    {
        if(bounds_cond[b].type_x == boundary_type::PRESSURE)
            pressure_cond<component::X>(mesh, b, bounds_cond[b].func_x, f);
        if(bounds_cond[b].type_y == boundary_type::PRESSURE)
            pressure_cond<component::Y>(mesh, b, bounds_cond[b].func_y, f);
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(bounds_cond[b].type_x == boundary_type::DISPLACEMENT)
            displacement_cond<component::X, MatrixMajor>(mesh, b, temperature_nodes, bounds_cond[b].func_x, K_bound, f);
        if(bounds_cond[b].type_y == boundary_type::DISPLACEMENT)
            displacement_cond<component::Y, MatrixMajor>(mesh, b, temperature_nodes, bounds_cond[b].func_y, K_bound, f);
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(bounds_cond[b].type_x == boundary_type::DISPLACEMENT)
            for(auto node : temperature_nodes[b])
                f[2*node]   = bounds_cond[b].func_x(mesh.coord(node, 0), mesh.coord(node, 1));

        if(bounds_cond[b].type_y == boundary_type::DISPLACEMENT)
            for(auto node : temperature_nodes[b])
                f[2*node+1] = bounds_cond[b].func_y(mesh.coord(node, 0), mesh.coord(node, 1));
    }
}

template<class Type, class Index>
static void integrate_right_part(const mesh_2d<Type, Index> &mesh,
                                 const distributed_load<Type>& right_part,
                                 Eigen::Matrix<Type, Eigen::Dynamic, 1> &f)
{
    matrix<Type> coords, jacobi_matrices;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        approx_quad_nodes_coords(mesh, e, el, coords);
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            f[2*mesh.node_number(el, i)]   += integrate_right_part_function(e, i, coords, jacobi_matrices, right_part.func_x);
            f[2*mesh.node_number(el, i)+1] += integrate_right_part_function(e, i, coords, jacobi_matrices, right_part.func_y);
        }
    }
}

template<class Type, class Index>
static Eigen::Matrix<Type, Eigen::Dynamic, 1> 
    stationary(const mesh_2d<Type, Index> &mesh, const parameters<Type> &params,
               const std::vector<boundary_condition<Type>> &bounds_cond,
               const distributed_load<Type>& right_part,
               const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(2*mesh.nodes_count());
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                                      K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    
    double time = omp_get_wtime();
    create_matrix(mesh, params, bounds_cond, K, K_bound, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    integrate_right_part(mesh, right_part, f);

    time = omp_get_wtime();
    boundary_condition_calc(mesh, kinematic_nodes_vectors(mesh, bounds_cond), bounds_cond, K_bound, f);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;

    return u;
}

template<class Type, class Index>
static std::array<std::vector<Type>, 6>
    strains_and_stress_loc(const mesh_2d<Type, Index>& mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1>& u, const std::array<Type, 3>& D)
{
    std::vector<Type> eps11  (mesh.nodes_count(), 0.0),
                      eps22  (mesh.nodes_count(), 0.0),
                      eps12  (mesh.nodes_count(), 0.0),
                      sigma11(mesh.nodes_count(), 0.0),
                      sigma22(mesh.nodes_count(), 0.0),
                      sigma12(mesh.nodes_count(), 0.0);

    eps11.shrink_to_fit();

    std::array<Type, 4> jacobi;
    std::array<Type, 3> loc_eps;
    std::vector<uint8_t> repeating(mesh.nodes_count(), 0);
    const finite_element::element_2d_integrate_base<Type> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            ++repeating[mesh.node_number(el, i)];
            const std::array<Type, 2>& node = e->node(i);
            memset(jacobi.data(), 0, jacobi.size() * sizeof(Type));
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                jacobi[0] += mesh.coord(mesh.node_number(el, j), 0) * e->Nxi (j, node[0], node[1]);
                jacobi[1] += mesh.coord(mesh.node_number(el, j), 0) * e->Neta(j, node[0], node[1]);
                jacobi[2] += mesh.coord(mesh.node_number(el, j), 1) * e->Nxi (j, node[0], node[1]);
                jacobi[3] += mesh.coord(mesh.node_number(el, j), 1) * e->Neta(j, node[0], node[1]);
            }

            memset(loc_eps.data(), 0, loc_eps.size() * sizeof(Type));
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                const Type jacobian = jacobi[0]*jacobi[3] - jacobi[1]*jacobi[2],
                           dx1 =  jacobi[3] * e->Nxi(j, node[0], node[1]) - jacobi[2] * e->Neta(j, node[0], node[1]),
                           dx2 = -jacobi[1] * e->Nxi(j, node[0], node[1]) + jacobi[0] * e->Neta(j, node[0], node[1]);
                loc_eps[0] +=  dx1 * u[2*mesh.node_number(el, j)  ]  / jacobian;
                loc_eps[1] +=  dx2 * u[2*mesh.node_number(el, j)+1]  / jacobian;
                loc_eps[2] += (dx2 * u[2*mesh.node_number(el, j)  ] +
                               dx1 * u[2*mesh.node_number(el, j)+1]) / jacobian;
            }

            eps11  [mesh.node_number(el, i)] += loc_eps[0];
            eps22  [mesh.node_number(el, i)] += loc_eps[1];
            eps12  [mesh.node_number(el, i)] += loc_eps[2];
            sigma11[mesh.node_number(el, i)] += D[0] * loc_eps[0] + D[1] * loc_eps[1];
            sigma22[mesh.node_number(el, i)] += D[1] * loc_eps[0] + D[0] * loc_eps[1];
            sigma12[mesh.node_number(el, i)] += D[2] * loc_eps[2];
        }
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        eps11  [i] /=   repeating[i];
        eps22  [i] /=   repeating[i];
        eps12  [i] /= 2*repeating[i];
        sigma11[i] /=   repeating[i];
        sigma22[i] /=   repeating[i];
        sigma12[i] /= 2*repeating[i];
    }

    return {std::move(eps11), std::move(eps22), std::move(eps12), std::move(sigma11), std::move(sigma22), std::move(sigma12)};
}

template<class Type, class Index>
static std::array<std::vector<Type>, 3>
    approx_all_eps_in_all_quad(const mesh_2d<Type, Index> &mesh, const std::vector<Index> &shifts,
                               const std::vector<Type> &eps11, const std::vector<Type> &eps22, const std::vector<Type> &eps12)
{
    std::vector<Type> all_eps11(shifts.back(), 0.), all_eps22(shifts.back(), 0.), all_eps12(shifts.back(), 0.);
    const finite_element::element_2d_integrate_base<Type> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t q = 0, shift = shifts[el]; q < e->qnodes_count(); ++q, ++shift)
            for(size_t i = 0; i < e->nodes_count(); ++i)
            {
                all_eps11[shift] += eps11[mesh.node_number(el, i)] * e->qN(i, q);
                all_eps22[shift] += eps22[mesh.node_number(el, i)] * e->qN(i, q);
                all_eps12[shift] += eps12[mesh.node_number(el, i)] * e->qN(i, q);
            }
    }
    return {std::move(all_eps11), std::move(all_eps22), std::move(all_eps12)};
}

template<class Type, class Index>
static void stress_nonloc(const mesh_2d<Type, Index> &mesh, const std::array<Type, 3> &D,
                          const std::vector<Type> &eps11,   const std::vector<Type> &eps22,   const std::vector<Type> &eps12,
                                std::vector<Type> &sigma11,       std::vector<Type> &sigma22,       std::vector<Type> &sigma12,
                          const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    const Type p2 = 1. - p1;
    const finite_element::element_2d_integrate_base<Type> *eNL = nullptr;
    const std::vector<Index> shifts_quad = quadrature_shifts_init(mesh);
    const matrix<Type> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts_quad);
    const matrix<Type> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);
    auto [all_eps11, all_eps22, all_eps12] = approx_all_eps_in_all_quad(mesh, shifts_quad, eps11, eps22, eps12);
    for(size_t node = 0; node < mesh.nodes_count(); ++node)
        for(const auto elNL : mesh.neighbor(node))
        {
            eNL = mesh.element_2d(mesh.element_type(elNL));
            for(size_t q = 0, shift = shifts_quad[elNL]; q < eNL->qnodes_count(); ++q, ++shift)
            {
                const Type finit = influence_fun(mesh.coord(node, 0), all_quad_coords(shift, 0), mesh.coord(node, 1), all_quad_coords(shift, 1)) *
                                   (all_jacobi_matrices(shift, 0)*all_jacobi_matrices(shift, 3) - all_jacobi_matrices(shift, 1)*all_jacobi_matrices(shift, 2));
                sigma11[node] += p2 * finit * (D[0] * all_eps11[shift] + D[1] * all_eps22[shift]);
                sigma22[node] += p2 * finit * (D[1] * all_eps11[shift] + D[0] * all_eps22[shift]);
                sigma12[node] += p2 * finit *  D[2] * all_eps12[shift];
            }
        }
}

template<class Type, class Index>
static std::array<std::vector<Type>, 6>
    strains_and_stress(const mesh_2d<Type, Index> &mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1> &u, const parameters<Type> &params,
                       const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    static constexpr Type MAX_LOCAL_WEIGHT = 0.999;
    bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
    const std::array<Type, 3> D = hooke_matrix(params);
    auto [eps11, eps22, eps12, sigma11, sigma22, sigma12] = strains_and_stress_loc(mesh, u, D);

    if(nonlocal)
    {
        for(size_t i = 0; i < mesh.nodes_count(); ++i)
        {
            sigma11[i] *= p1;
            sigma22[i] *= p1;
            sigma12[i] *= p1;
        }
        stress_nonloc(mesh, D, eps11, eps22, eps12, sigma11, sigma22, sigma12, p1, influence_fun);
    }

    return {std::move(eps11), std::move(eps22), std::move(eps12), std::move(sigma11), std::move(sigma22), std::move(sigma12)};
}

template<class Type, class Index>
void raw_output(const std::string &path,          const mesh_2d<Type, Index> &mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1> &u,
                const std::vector<Type> &eps11,   const std::vector<Type> &eps22,   const std::vector<Type> &eps12,
                const std::vector<Type> &sigma11, const std::vector<Type> &sigma22, const std::vector<Type> &sigma12)
{
    std::ofstream fout_ux(path + std::string("u_x.csv")),
                  fout_uy(path + std::string("u_y.csv"));
    fout_ux.precision(20);
    fout_uy.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        fout_ux << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i) << std::endl;
        fout_uy << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i+1) << std::endl;
    }

    std::ofstream fout_eps11(path + std::string("eps11.csv")),
                  fout_eps22(path + std::string("eps22.csv")),
                  fout_eps12(path + std::string("eps12.csv"));
    fout_eps11.precision(20);
    fout_eps22.precision(20);
    fout_eps12.precision(20);
    for(size_t i = 0; i < eps11.size(); ++i)
        fout_eps11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps11[i] << std::endl;

    for(size_t i = 0; i < eps22.size(); ++i)
        fout_eps22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps22[i] << std::endl;

    for(size_t i = 0; i < eps12.size(); ++i)
        fout_eps12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps12[i] << std::endl;

    std::ofstream fout_sigma11(path + std::string("sigma11.csv")),
                  fout_sigma22(path + std::string("sigma22.csv")),
                  fout_sigma12(path + std::string("sigma12.csv"));
    fout_sigma11.precision(20);
    fout_sigma22.precision(20);
    fout_sigma12.precision(20);
    for(size_t i = 0; i < sigma11.size(); ++i)
        fout_sigma11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma11[i] << std::endl;

    for(size_t i = 0; i < sigma22.size(); ++i)
        fout_sigma22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma22[i] << std::endl;

    for(size_t i = 0; i < sigma12.size(); ++i)
        fout_sigma12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma12[i] << std::endl;
}

template<class Type, class Index>
void save_as_vtk(const std::string &path,          const mesh_2d<Type, Index> &mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1> &u,
                 const std::vector<Type> &eps11,   const std::vector<Type> &eps22,   const std::vector<Type> &eps12,
                 const std::vector<Type> &sigma11, const std::vector<Type> &sigma22, const std::vector<Type> &sigma12)
{
    std::ofstream fout(path);
    fout.precision(20);

    fout << "# vtk DataFile Version 4.2" << std::endl
         << "Temperature"                << std::endl
         << "ASCII"                      << std::endl
         << "DATASET UNSTRUCTURED_GRID"  << std::endl;

    fout << "POINTS " << mesh.nodes_count() << " double" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << mesh.coord(i, 0) << " " << mesh.coord(i, 1) << " 0" << std::endl;

    fout << "CELLS " << mesh.elements_count() << " " << mesh.elements_count() * (mesh.element_type(0) == 3 ? 5 : 9) << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
    {
        if (mesh.element_type(i) == 3)
        {
            fout << 4 << " " << mesh.node_number(i, 0) << " "
                             << mesh.node_number(i, 1) << " "
                             << mesh.node_number(i, 2) << " "
                             << mesh.node_number(i, 3);
        }
        else
        {
            fout << 8 << " " << mesh.node_number(i, 0) << " "
                             << mesh.node_number(i, 2) << " "
                             << mesh.node_number(i, 4) << " "
                             << mesh.node_number(i, 6) << " "
                             << mesh.node_number(i, 1) << " "
                             << mesh.node_number(i, 3) << " "
                             << mesh.node_number(i, 5) << " "
                             << mesh.node_number(i, 7);
        }
        fout << std::endl;
    }
        //fout << 4 << " " << mesh.node_number(i, 0) << " "
        //                 << mesh.node_number(i, 1) << " "
        //                 << mesh.node_number(i, 2) << " "
        //                 << mesh.node_number(i, 3) << std::endl;

    fout << "CELL_TYPES " << mesh.elements_count() << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        fout << (mesh.element_type(i) == 3 ? 9 : 23) << std::endl;

    fout << "POINT_DATA " << mesh.nodes_count() << std::endl;

    fout << "SCALARS U_X double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << u[2*i] << std::endl;

    fout << "SCALARS U_Y double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << u[2*i+1] << std::endl;

    fout << "SCALARS EPS_XX double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps11.size(); ++i)
        fout << eps11[i] << std::endl;

    fout << "SCALARS EPS_YY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps22.size(); ++i)
        fout << eps22[i] << std::endl;

    fout << "SCALARS EPS_XY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < eps12.size(); ++i)
        fout << eps12[i] << std::endl;

    fout << "SCALARS SIGMA_XX double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma11.size(); ++i)
        fout << sigma11[i] << std::endl;

    fout << "SCALARS SIGMA_YY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma22.size(); ++i)
        fout << sigma22[i] << std::endl;

    fout << "SCALARS SIGMA_XY double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < sigma12.size(); ++i)
        fout << sigma12[i] << std::endl;
}

}

#endif