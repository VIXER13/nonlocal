#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <tuple>
#include <algorithm>
#include "omp.h"
#include "finite_element_routine.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/PardisoSupport"

namespace heat_equation_with_nonloc
{

enum class boundary_type : uint8_t {TEMPERATURE, FLOW};

template<class Type>
struct boundary_condition
{
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    std::function<Type(Type, Type)> func = [](Type, Type) { return 0.; };
    boundary_type type = boundary_type::FLOW;
};

template<class Type>
static Type integrate_basic(const metamath::finite_element::element_2d_integrate_base<Type> *const e,
                            const size_t i, const matrix<Type>& jacobi_matrices)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->nodes_count(); ++q)
        integral += e->weight(q) * e->qN(i, q) *
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2));
    return integral;
}

template<class Type>
static Type integrate_basic_pair(const metamath::finite_element::element_2d_integrate_base<Type> *const e,
                                 const size_t i, const size_t j, const matrix<Type>& jacobi_matrices, size_t shift = 0)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->nodes_count(); ++q, ++shift)
        integral += e->weight(q) * e->qN(i, q) * e->qN(j, q) *
                    (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2));
    return integral;
}

template<class Type>
static Type integrate_gradient_pair(const metamath::finite_element::element_2d_integrate_base<Type> *const e,
                                    const size_t i, const size_t j, const matrix<Type>& jacobi_matrices, size_t shift = 0)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
        integral += (( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                     ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2)) +
                     (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                     (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0))) /
                    (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
    return integral;
}

// Influence_Function - функтор с сигнатурой Type(Type, Type, Type, Type)
template<class Type, class Finite_Element_2D_Pointer, class Influence_Function>
static Type integrate_gradient_pair_nonloc(const Finite_Element_2D_Pointer& eL, const Finite_Element_2D_Pointer& eNL,
                                           const size_t iL, const size_t jNL, size_t shiftL, size_t shiftNL,
                                           const matrix<Type>& coords, const matrix<Type>& jacobi_matrices,
                                           const Influence_Function& influence_fun)
{
    const size_t sub = shiftNL;
    Type integral = 0, int_with_weight_x = 0, int_with_weight_y = 0, finit = 0;
    for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL)
    {
        int_with_weight_x = 0;
        int_with_weight_y = 0;
        for(size_t qNL = 0, shiftNL = sub; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL)
        {
            finit = eNL->weight(qNL) * influence_fun(coords(shiftL, 0), coords(shiftNL, 0), coords(shiftL, 1), coords(shiftNL, 1));
            int_with_weight_x += finit * ( eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 3) - eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 2));
            int_with_weight_y += finit * (-eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 1) + eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 0));
        }
        integral += eL->weight(qL) *
                    (int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)) +
                     int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)));
    }
    return integral;
}

template<class Type, class Index>
Type integrate_solution(const mesh_2d<Type, Index>& mesh, const Eigen::Matrix<Type, Eigen::Dynamic, 1>& x)
{
    Type integral = 0.;
    matrix<Type> jacobi_matrices;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                integral += e->weight(q) * e->qN(i, q) * x[mesh.node_number(el, i)] *
                            (jacobi_matrices(q, 0) * jacobi_matrices(q, 3) - jacobi_matrices(q, 1) * jacobi_matrices(q, 2));
    }
    return integral;
}

template<class Type, class Index, class Right_Part>
static void integrate_right_part(const mesh_2d<Type, Index>& mesh, const Right_Part& right_part, Eigen::Matrix<Type, Eigen::Dynamic, 1>& f)
{
    matrix<Type> coords, jacobi_matrices;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        approx_quad_nodes_coords(mesh, e, el, coords);
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            f[mesh.node_number(el, i)] += integrate_right_part_function(e, i, coords, jacobi_matrices, right_part);
    }
}

template<class Type, class Index>
static std::vector<bool> inner_nodes_vector(const mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond)
{
    std::vector<bool> inner_nodes(mesh.nodes_count(), true);
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(bounds_cond[b].type == boundary_type::TEMPERATURE)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                inner_nodes[*node] = false;
    return std::move(inner_nodes);
}

// Создаёт массив векторов, в каждом из которых хранятся номера узлов в которых заданы граничные условия первого рода на соответствующей границе,
// с тем условием, что среди всех векторов нет повторяющихся номеров. То есть номер узла на стыке будет записан лишь в один из векторов.
template<class Type, class Index>
static std::vector<std::vector<Index>>
    temperature_nodes_vectors(const mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond)
{
    std::vector<std::vector<Index>> temperature_nodes(bounds_cond.size());
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(bounds_cond[b].type == boundary_type::TEMPERATURE)
            for(auto [node, k] = std::make_tuple(mesh.boundary(b).cbegin(), size_t(0)); node != mesh.boundary(b).cend(); ++node)
            {
                for(k = 0; k < temperature_nodes.size(); ++k)
                    if(std::find(temperature_nodes[k].cbegin(), temperature_nodes[k].cend(), *node) != temperature_nodes[k].cend())
                        break;
                if(k == temperature_nodes.size())
                   temperature_nodes[b].push_back(*node);
            }
    return std::move(temperature_nodes);
}

template<class Type, class Index>
static std::array<std::vector<Index>, 4>
    mesh_analysis(const mesh_2d<Type, Index>& mesh, const std::vector<bool>& inner_nodes, const bool nonlocal)
{
    std::vector<Index> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0),
                       shifts_nonloc, shifts_bound_nonloc;

    mesh_run_loc(mesh,
        [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc](size_t i, size_t j, size_t el)
        { 
            if(mesh.node_number(el, i) >= mesh.node_number(el, j))
            {
                if(inner_nodes[mesh.node_number(el, i)] && inner_nodes[mesh.node_number(el, j)])
                    ++shifts_loc[el+1];
                else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                    ++shifts_bound_loc[el+1];
            }
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

        mesh_run_nonloc(mesh, 
            [&mesh, &inner_nodes, &shifts_nonloc, &shifts_bound_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
            { 
                if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL))
                {
                    if(inner_nodes[mesh.node_number(elL, iL)] && inner_nodes[mesh.node_number(elNL, jNL)])
                        ++shifts_nonloc[elL+1];
                    else if(mesh.node_number(elL, iL) != mesh.node_number(elNL, jNL))
                        ++shifts_bound_nonloc[elL+1];
                }
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

// Integrate_Rule - функтор с сигнатурой Type(const metamath::finite_element::element_2d_integrate_base<Type>*, 
//                                            const size_t, const size_t, const matrix<Type>&, size_t)
// Influence_Function - функтор с сигнатурой Type(Type, Type, Type, Type)
template<class Type, class Index, class Integrate_Rule, class Influence_Function>
static std::array<std::vector<Eigen::Triplet<Type, Index>>, 2>
    triplets_fill(const mesh_2d<Type, Index>& mesh, const std::vector<bool>& inner_nodes, const Integrate_Rule& integrate_rule,
                  const Type p1, const Influence_Function& influence_fun)
{
    static constexpr Type MAX_LOCAL_WEIGHT = 0.999;
    bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
    auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, inner_nodes, nonlocal);
    std::vector<Eigen::Triplet<Type, Index>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                             triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
    for(size_t i = 0, j = 0; i < inner_nodes.size(); ++i)
        if(!inner_nodes[i])
            triplets[j++] = Eigen::Triplet<Type, Index>(i, i, 1.);

    const std::vector<int> shifts_quad = quadrature_shifts_init(mesh);
    const matrix<Type> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);

    mesh_run_loc(mesh,
        [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_loc, &shifts_bound_loc,
         &shifts_quad, &all_jacobi_matrices, &integrate_rule, p1] (size_t i, size_t j, size_t el)
        {
            if(mesh.node_number(el, i) >= mesh.node_number(el, j))
            {
                const Type integral = p1 * integrate_rule(mesh.element_2d(mesh.element_type(el)), i, j, all_jacobi_matrices, shifts_quad[el]);
                if(inner_nodes[mesh.node_number(el, i)] && inner_nodes[mesh.node_number(el, j)])
                    triplets[shifts_loc[el]++] = Eigen::Triplet<Type>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
                else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                    triplets_bound[shifts_bound_loc[el]++] = inner_nodes[mesh.node_number(el, j)] ?
                                                             Eigen::Triplet<Type, Index>(mesh.node_number(el, j), mesh.node_number(el, i), integral) :
                                                             Eigen::Triplet<Type, Index>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
            }
        });

    if(nonlocal)
    {
        const matrix<Type> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts_quad);
        mesh_run_nonloc(mesh, 
            [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
             &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1]
            (size_t iL, size_t jNL, size_t elL, size_t elNL)
            {
                if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL))
                {
                    const Type integral = p2 * integrate_gradient_pair_nonloc(mesh.element_2d(mesh.element_type(elL )),
                                                                              mesh.element_2d(mesh.element_type(elNL)), 
                                                                              iL, jNL, shifts_quad[elL], shifts_quad[elNL],
                                                                              all_quad_coords, all_jacobi_matrices, influence_fun);
                    if(inner_nodes[mesh.node_number(elL, iL)] && inner_nodes[mesh.node_number(elNL, jNL)])
                        triplets[shifts_nonloc[elL]++] = Eigen::Triplet<Type>(mesh.node_number(elL, iL), mesh.node_number(elNL, jNL), integral);
                    else if(mesh.node_number(elL, iL) != mesh.node_number(elNL, jNL))
                        triplets_bound[shifts_bound_nonloc[elL]++] = inner_nodes[mesh.node_number(elNL, jNL)] ?
                                                                     Eigen::Triplet<Type, Index>(mesh.node_number(elNL, jNL), mesh.node_number(elL,   iL), integral) :
                                                                     Eigen::Triplet<Type, Index>(mesh.node_number(elL,  iL ), mesh.node_number(elNL, jNL), integral);
                }
            });
    }

    return {std::move(triplets), std::move(triplets_bound)};
}

// Вычисление матрицы теплопроводности (или теплоёмкости, в зависимости от функции integrate_rule)
// с разделением матрицы на две части, одна из которых является матрицей после учёта гран условий первого рода,
// а вторая та часть матрицы, которая уйдёт в правую часть при учёте гран условий первого рода.
// Integrate_Rule - функтор с сигнатурой Type(const metamath::finite_element::element_2d_integrate_base<Type>*, 
//                                            const size_t, const size_t, const matrix<Type>&, size_t)
template<int MatrixMajor, class Type, class Index, class Integrate_Rule>
static void create_matrix(const mesh_2d<Type, Index>& mesh,
                          const std::vector<boundary_condition<Type>>& bounds_cond,
                          const Integrate_Rule& integrate_rule, const Type p1, const std::function<Type(Type, Type, Type, Type)>& influence_fun,
                          Eigen::SparseMatrix<Type, MatrixMajor, Index>& K, Eigen::SparseMatrix<Type, MatrixMajor, Index>& K_bound)
{
    double time = omp_get_wtime();
    auto [triplets, triplets_bound] = triplets_fill(mesh, inner_nodes_vector(mesh, bounds_cond), integrate_rule, p1, influence_fun);
    std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;
    std::cout << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;

    K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
    triplets_bound.reserve(0);
    K.setFromTriplets(triplets.cbegin(), triplets.cend());
    std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
}

template<class Type, int MatrixMajor, class Index>
static void boundary_condition_calc(const mesh_2d<Type, Index>& mesh, const std::vector<std::vector<Index>>& temperature_nodes,
                                    const std::vector<boundary_condition<Type>>& bounds_cond,
                                    const Type tau, const Eigen::SparseMatrix<Type, MatrixMajor, Index>& K_bound,
                                    Eigen::Matrix<Type, Eigen::Dynamic, 1>& f)
{
    matrix<Type> coords, jacobi_matrices;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(bounds_cond[b].type == boundary_type::FLOW)
            for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
            {
                const auto& be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
                approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
                approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
                for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
                    f[mesh.boundary(b)(el, i)] += tau*integrate_boundary_gradient(be, i, coords, jacobi_matrices, bounds_cond[b].func);
            }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(const auto node : temperature_nodes[b])
        {
            const Type temp = bounds_cond[b].func(mesh.coord(node, 0), mesh.coord(node, 1));
            for(typename Eigen::SparseMatrix<Type>::InnerIterator it(K_bound, node); it; ++it)
                f[it.row()] -= temp * it.value();
        }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(const auto node : temperature_nodes[b])
            f[node] = bounds_cond[b].func(mesh.coord(node, 0), mesh.coord(node, 1));
}

// Нелокальное условие для задачи Неймана.
template<int MatrixMajor, class Type, class Index>
static Eigen::SparseMatrix<Type, MatrixMajor, Index> nonlocal_condition(const mesh_2d<Type, Index>& mesh)
{
    size_t triplets_count = 0;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
        triplets_count += mesh.element_2d(mesh.element_type(el))->nodes_count();

    matrix<Type> jacobi_matrices;
    std::vector<Eigen::Triplet<Type, Index>> triplets(triplets_count);

    triplets_count = 0;
    const metamath::finite_element::element_2d_integrate_base<Type> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            triplets[triplets_count++] = Eigen::Triplet<Type, Index>(mesh.nodes_count(), mesh.node_number(el, i), 
                                                                     integrate_basic(e, i, jacobi_matrices));
    }

    Eigen::SparseMatrix<Type, MatrixMajor, Index> K_last_row(mesh.nodes_count()+1, mesh.nodes_count()+1);
    K_last_row.setFromTriplets(triplets.cbegin(), triplets.cend());
    return K_last_row;
}

// Решает стационарное уравнение теплопроводности
template<class Type, class Index>
Eigen::Matrix<Type, Eigen::Dynamic, 1>
    stationary(const mesh_2d<Type, Index>& mesh,
               const std::vector<boundary_condition<Type>>& bounds_cond,
               const std::function<Type(Type, Type)>& right_part,
               const Type p1, const std::function<Type(Type, Type, Type, Type)>& influence_fun,
               const Type volume)
{
    const bool neumann_task = std::all_of(bounds_cond.cbegin(), bounds_cond.cend(),
                              [](const boundary_condition<Type> &bound) { return bound.type == boundary_type::FLOW; });

    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count(),
                                                              neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count()),
                                                      K_bound(neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count(),
                                                              neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count());
    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count());
    if(neumann_task)
        f[mesh.nodes_count()] = volume;

    double time = omp_get_wtime();
    integrate_right_part(mesh, right_part, f);
    std::cout << "Right part Integrate: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    create_matrix(mesh, bounds_cond, integrate_gradient_pair<Type>, p1, influence_fun, K, K_bound);
    if(neumann_task)
        K += nonlocal_condition<Eigen::ColMajor>(mesh);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;
    std::cout << "Nonzero elements count = " << K.nonZeros() + K_bound.nonZeros() << std::endl;

    time = omp_get_wtime();
    boundary_condition_calc(mesh, temperature_nodes_vectors(mesh, bounds_cond), bounds_cond, Type(1.), K_bound, f);
    std::cout << "Boundary filling: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> T;
    if(neumann_task)
    {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
        solver.compute(K);
        T = solver.solve(f);
    }
    else
    {
        Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
        solver.compute(K);
        T = solver.solve(f);
    }
    
    std::cout << "System solving: " << omp_get_wtime() - time << std::endl;
    return T;
}

template<class Type, class Index>
void nonstationary(const std::string& path,
                   const mesh_2d<Type, Index>& mesh, const Type tau, const size_t time_steps,
                   const std::vector<boundary_condition<Type>>& bounds_cond,
                   const std::function<Type(Type, Type)>& init_dist,
                   const std::function<Type(Type, Type)>& right_part,
                   const Type p1, const std::function<Type(Type, Type, Type, Type)>& influence_fun,
                   const uint64_t print_frequency)
{
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (mesh.nodes_count(), mesh.nodes_count()),
                                                      K_bound(mesh.nodes_count(), mesh.nodes_count()),
                                                      C      (mesh.nodes_count(), mesh.nodes_count()),
                                                      C_bound(mesh.nodes_count(), mesh.nodes_count());
    create_matrix(mesh, bounds_cond, integrate_gradient_pair<Type>, p1, influence_fun, K, K_bound);
    create_matrix(mesh, bounds_cond, integrate_basic_pair<Type>,    1., influence_fun, C, C_bound);
    C_bound.setZero();
    K *= tau;
    K_bound *= tau;
    K += C;

    std::vector<std::vector<Index>> temperature_nodes = temperature_nodes_vectors(mesh, bounds_cond);
    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(auto node : temperature_nodes[b])
            K.coeffRef(node, node) = 1.;

    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(mesh.nodes_count()),
                                           T_prev(mesh.nodes_count()), T(mesh.nodes_count());
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        T_prev[i] = init_dist(mesh.coord(i, 0), mesh.coord(i, 1));

    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);

    if(print_frequency != uint64_t(-1))
    {
        mesh.print_to_file(path+std::string("0.csv"), T_prev.data());
        std::cout << "step = " << 0 << " Volume = " << integrate_solution(mesh, T_prev) << std::endl;
    }
    for(size_t i = 1; i < time_steps; ++i)
    {
        f = C.template selfadjointView<Eigen::Lower>() * T_prev;
        integrate_right_part(mesh, right_part, f);
        boundary_condition_calc(mesh, temperature_nodes, bounds_cond, tau, K_bound, f);
        T = solver.solve(f);
        if(i % print_frequency == 0)
        {
            mesh.print_to_file(path+std::to_string(i)+std::string(".csv"), T.data());
            std::cout << "step = " << i << " Volume = " << integrate_solution(mesh, T) << std::endl;
        }
        T_prev.swap(T);
    }
}

template<class Type, class Index, class Vector>
static void raw_output(const std::string& path, const mesh_2d<Type, Index>& mesh, const Vector& T)
{
    std::ofstream fout(path + std::string{"T.csv"});
    fout.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << T[i] << std::endl;
}

template<class Type, class Index, class Vector>
static void save_as_vtk(const std::string& path, const mesh_2d<Type, Index>& mesh, const Vector& T)
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

    fout << "CELLS " << mesh.elements_count() << " " << mesh.elements_count() * 5 << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        fout << 4 << " " << mesh.node_number(i, 0) << " "
                         << mesh.node_number(i, 1) << " "
                         << mesh.node_number(i, 2) << " "
                         << mesh.node_number(i, 3) << std::endl;

    fout << "CELL_TYPES " << mesh.elements_count() << std::endl;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        fout << 9 << std::endl;

    fout << "POINT_DATA " << mesh.nodes_count() << std::endl;

    fout << "SCALARS TEMPERATURE double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        fout << T[i] << std::endl;
}

}

#endif