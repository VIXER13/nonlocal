#include <set>
#include <algorithm>
#include "omp.h"
#include "Eigen/Sparse"
#include "finite_element_routine.hpp"
#include "heat_equation_solver.hpp"

namespace heat_equation_with_nonloc
{

[[maybe_unused]] static void save_as_vtk(const std::string &path, const mesh_2d<double> &mesh, const Eigen::VectorXd &sol)
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

    for(size_t i = 0; i < size_t(sol.size()); ++i)
        fout << sol[i] << std::endl;
}

template<class Type>
static Type integrate_basic(const finite_element::element_2d_integrate_base<Type> *const e,
                            const size_t i, const matrix<Type> &jacobi_matrices)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->nodes_count(); ++q)
        integral += e->weight(q) * e->qN(i, q) *
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2));
    return integral;
}

template<class Type>
static Type integrate_basic_pair(const finite_element::element_2d_integrate_base<Type> *const e,
                                 const size_t i, const size_t j, const matrix<Type> &jacobi_matrices, size_t shift = 0)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->nodes_count(); ++q, ++shift)
        integral += e->weight(q) * e->qN(i, q) * e->qN(j, q) *
                    (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2));
    return integral;
}

template<class Type>
static Type integrate_gradient_pair(const finite_element::element_2d_integrate_base<Type> *const e,
                                    const size_t i, const size_t j, const matrix<Type> &jacobi_matrices, size_t shift = 0)
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

template<class Type>
static Type integrate_gradient_pair_nonloc(const finite_element::element_2d_integrate_base<Type> *const eL,
                                           const finite_element::element_2d_integrate_base<Type> *const eNL,
                                           const size_t iL, const size_t jNL, size_t shiftL, size_t shiftNL,
                                           const matrix<Type> &coords, const matrix<Type> &jacobi_matrices,
                                           const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    const size_t sub = shiftNL;
    Type integral = 0., int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
    for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL)
    {
        int_with_weight_x = 0.;
        int_with_weight_y = 0.;
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

template<class Type>
static Type integrate_function(const finite_element::element_2d_integrate_base<Type> *const e, const size_t i,
                               const matrix<Type> &coords, const matrix<Type> &jacobi_matrices,
                               const std::function<Type(Type, Type)> &fun)
{
    Type integral = 0.;
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        integral += e->weight(q) * e->qN(i, q) * fun(coords(q, 0), coords(q, 1)) *
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2));
    return integral;
}

template<class Type>
static Type integrate_flow_bound(const finite_element::element_1d_integrate_base<Type> *const be, const size_t i,
                                 const matrix<Type> &coords, const matrix<Type> &jacobi_matrices, 
                                 const std::function<Type(Type, Type)> &fun)
{
    Type integral = 0.;
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        integral += fun(coords(q, 0), coords(q, 1)) * be->weight(q) * be->qN(i, q) *
                    sqrt(jacobi_matrices(q, 0)*jacobi_matrices(q, 0) + jacobi_matrices(q, 1)*jacobi_matrices(q, 1));
    return integral;
}

double integrate_solution(const mesh_2d<double> &mesh, const Eigen::VectorXd &x)
{
    double integral = 0.;
    matrix<double> jacobi_matrices;
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                integral += e->weight(q) * e->qN(i, q) * x[mesh.node_number(el, i)] *
                            (jacobi_matrices(q, 0) * jacobi_matrices(q, 3) - jacobi_matrices(q, 1) * jacobi_matrices(q, 2));
    }
    return integral;
}

// Создаёт std::set из номеров узлов, на которых заданы граничные условия первого рода.
// Данный std::set необходим для заполнения информации о сетке в функции mesh_analysis.
static std::set<uint32_t> temperature_nodes_set(const mesh_2d<double> &mesh,
                              const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::set<uint32_t> temperature_nodes;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(std::get<boundary_type>(bounds_cond[b]) == boundary_type::TEMPERATURE)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                temperature_nodes.insert(*node);
    return temperature_nodes;
}

// Создаёт массив векторов, в каждом из которых хранятся номера узлов в которых заданы граничные условия первого рода на соответствующей границе,
// с тем условием, что среди всех векторов нет повторяющихся номеров. То есть номер узла на стыке будет записан лишь в один из векторов.
static std::vector<std::vector<uint32_t>> temperature_nodes_vectors(const mesh_2d<double> &mesh,
                                              const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::vector<std::vector<uint32_t>> temperature_nodes(bounds_cond.size());
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(std::get<boundary_type>(bounds_cond[b]) == boundary_type::TEMPERATURE)
            for(auto [node, k] = std::make_tuple(mesh.boundary(b).cbegin(), size_t(0)); node != mesh.boundary(b).cend(); ++node)
            {
                for(k = 0; k < temperature_nodes.size(); ++k)
                    if(std::find(temperature_nodes[k].cbegin(), temperature_nodes[k].cend(), *node) != temperature_nodes[k].cend())
                        break;
                if(k == temperature_nodes.size())
                   temperature_nodes[b].push_back(*node);
            }
    return temperature_nodes;
}

static void temperature(const mesh_2d<double> &mesh, const Eigen::SparseMatrix<double> &K_bound, Eigen::VectorXd &f,
                        const std::function<double(double, double)> &boundaryFun, const size_t node)
{
    const double temp = boundaryFun(mesh.coord(node, 0), mesh.coord(node, 1));
    for(typename Eigen::SparseMatrix<double>::InnerIterator it(K_bound, node); it; ++it)
        f[it.row()] -= temp * it.value();
}

static void boundary_condition(const mesh_2d<double> &mesh, const std::vector<std::vector<uint32_t>> &temperature_nodes,
                               const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                               const double tau, const Eigen::SparseMatrix<double> &K_bound, Eigen::VectorXd &f)
{
    const finite_element::element_1d_integrate_base<double> *be = nullptr;
    matrix<double> coords, jacobi_matrices;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(std::get<boundary_type>(bounds_cond[b]) == boundary_type::FLOW)
            for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
            {
                be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
                approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
                approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
                for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
                    f[mesh.boundary(b)(el, i)] += tau*integrate_flow_bound(be, i, coords, jacobi_matrices, std::get<1>(bounds_cond[b]));
            }


    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(auto node : temperature_nodes[b])
            temperature(mesh, K_bound, f, std::get<1>(bounds_cond[b]), node);

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(auto node : temperature_nodes[b])
            f[node] = std::get<1>(bounds_cond[b])(mesh.coord(node, 0), mesh.coord(node, 1));
}

static std::array<std::vector<uint32_t>, 4>
    mesh_analysis(const mesh_2d<double> &mesh, const std::set<uint32_t> &temperature_nodes, const bool nonlocal)
{
    std::vector<uint32_t> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0),
                          shifts_nonloc, shifts_bound_nonloc;

    mesh_run_loc(mesh,
        [&mesh, &temperature_nodes, &shifts_loc, &shifts_bound_loc](size_t i, size_t j, size_t el)
        { 
            if(mesh.node_number(el, i) >= mesh.node_number(el, j))
            {
                if(temperature_nodes.find(mesh.node_number(el, i)) == temperature_nodes.cend() &&
                   temperature_nodes.find(mesh.node_number(el, j)) == temperature_nodes.cend())
                    ++shifts_loc[el+1];
                else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                    ++shifts_bound_loc[el+1];
            }
        });

    shifts_loc[0] = temperature_nodes.size();
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
            [&mesh, &temperature_nodes, &shifts_nonloc, &shifts_bound_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
            { 
                if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL))
                {
                    if(temperature_nodes.find(mesh.node_number(elL , iL )) == temperature_nodes.cend() &&
                       temperature_nodes.find(mesh.node_number(elNL, jNL)) == temperature_nodes.cend())
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

static std::array<std::vector<Eigen::Triplet<double>>, 2>
    triplets_fill(const mesh_2d<double> &mesh, const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                  const double p1, const std::function<double(double, double, double, double)> &influence_fun,
                  const std::function<double(const finite_element::element_2d_integrate_base<double>*, 
                                             const size_t, const size_t, const matrix<double>&, size_t)> &integrate_rule)
{
    static constexpr double MAX_LOCAL_WEIGHT = 0.999;
    bool nonlocal = p1 < MAX_LOCAL_WEIGHT;

    const std::set<uint32_t> temperature_nodes = temperature_nodes_set(mesh, bounds_cond);
    auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, temperature_nodes, nonlocal);

    std::vector<Eigen::Triplet<double>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                        triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
    for(auto [it, i] = std::make_tuple(temperature_nodes.cbegin(), size_t(0)); it != temperature_nodes.cend(); ++it, ++i)
        triplets[i] = Eigen::Triplet<double>(*it, *it, 1.);

    const std::vector<uint32_t> shifts_quad = quadrature_shifts_init(mesh);
    const matrix<double> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);

    mesh_run_loc(mesh,
        [&mesh, &temperature_nodes, &triplets, &triplets_bound, &shifts_loc, &shifts_bound_loc,
         &shifts_quad, &all_jacobi_matrices, &integrate_rule, p1] (size_t i, size_t j, size_t el)
        {
            if(mesh.node_number(el, i) >= mesh.node_number(el, j))
            {
                const double integral = p1 * integrate_rule(mesh.element_2d(mesh.element_type(el)), i, j, all_jacobi_matrices, shifts_quad[el]);
                if(temperature_nodes.find(mesh.node_number(el, i)) == temperature_nodes.cend() &&
                   temperature_nodes.find(mesh.node_number(el, j)) == temperature_nodes.cend())
                    triplets[shifts_loc[el]++] = Eigen::Triplet<double>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
                else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                    triplets_bound[shifts_bound_loc[el]++] = temperature_nodes.find(mesh.node_number(el, j)) == temperature_nodes.cend() ?
                                                             Eigen::Triplet<double>(mesh.node_number(el, j), mesh.node_number(el, i), integral) :
                                                             Eigen::Triplet<double>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
            }
        });

    if(nonlocal)
    {
        const matrix<double> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts_quad);
        mesh_run_nonloc(mesh, 
            [&mesh, &temperature_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
             &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1]
            (size_t iL, size_t jNL, size_t elL, size_t elNL)
            {
                if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL))
                {
                    const double integral = p2 * integrate_gradient_pair_nonloc(mesh.element_2d(mesh.element_type(elL )),
                                                                                mesh.element_2d(mesh.element_type(elNL)), 
                                                                                iL, jNL, shifts_quad[elL], shifts_quad[elNL],
                                                                                all_quad_coords, all_jacobi_matrices, influence_fun);
                    if(temperature_nodes.find(mesh.node_number(elL , iL )) == temperature_nodes.cend() &&
                       temperature_nodes.find(mesh.node_number(elNL, jNL)) == temperature_nodes.cend())
                        triplets[shifts_nonloc[elL]++] = Eigen::Triplet<double>(mesh.node_number(elL, iL), mesh.node_number(elNL, jNL), integral);
                    else if(mesh.node_number(elL, iL) != mesh.node_number(elNL, jNL))
                        triplets_bound[shifts_bound_nonloc[elL]++] = temperature_nodes.find(mesh.node_number(elNL, jNL)) == temperature_nodes.cend() ?
                                                                     Eigen::Triplet<double>(mesh.node_number(elNL, jNL), mesh.node_number(elL,   iL), integral) :
                                                                     Eigen::Triplet<double>(mesh.node_number(elL,  iL ), mesh.node_number(elNL, jNL), integral);
                }
            });
    }

    return {std::move(triplets), std::move(triplets_bound)};
}

// Вычисление матрицы теплопроводности (или теплоёмкости, в зависимости от функции integrate_rule)
// с разделением матрицы на две части, одна из которых является матрицей после учёта гран условий первого рода,
// а вторая та часть матрицы, которая уйдёт в правую часть при учёте гран условий первого рода.
static void create_matrix(const mesh_2d<double> &mesh,
                          const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                          Eigen::SparseMatrix<double> &K, Eigen::SparseMatrix<double> &K_bound,
                          const std::function<double(const finite_element::element_2d_integrate_base<double>*, 
                                                     const size_t, const size_t, const matrix<double>&, size_t)> &integrate_rule,
                          const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    double time = omp_get_wtime();
    auto [triplets, triplets_bound] = triplets_fill(mesh, bounds_cond, p1, influence_fun, integrate_rule);
    std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;
    std::cout << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;

    K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
    triplets_bound.clear();
    K.setFromTriplets(triplets.cbegin(), triplets.cend());
    std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
}

static void integrate_right_part(const mesh_2d<double> &mesh, Eigen::VectorXd &f, const std::function<double(double, double)> &right_part)
{
    matrix<double> coords, jacobi_matrices;
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        approx_quad_nodes_coords(mesh, e, el, coords);
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            f[mesh.node_number(el, i)] += integrate_function(e, i, coords, jacobi_matrices, right_part);
    }
}

// Нелокальное условие для задачи Неймана.
static Eigen::SparseMatrix<double> nonlocal_condition(const mesh_2d<double> &mesh)
{
    size_t triplets_count = 0;
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ++triplets_count;
    }

    matrix<double> jacobi_matrices;
    std::vector<Eigen::Triplet<double>> triplets(triplets_count);

    triplets_count = 0;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            triplets[triplets_count++] = Eigen::Triplet<double>(mesh.nodes_count(), mesh.node_number(el, i), 
                                                                integrate_basic(e, i, jacobi_matrices));
    }

    Eigen::SparseMatrix<double> K_last_row(mesh.nodes_count()+1, mesh.nodes_count()+1);
    K_last_row.setFromTriplets(triplets.begin(), triplets.end());
    return K_last_row;
}

// Решает стационарное уравнение теплопроводности
void stationary(const std::string &path, const mesh_2d<double> &mesh,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                const std::function<double(double, double)> &right_part,
                const double p1, const std::function<double(double, double, double, double)> &influence_fun,
                const double volume)
{
    Eigen::SparseMatrix<double> K, K_bound;
    Eigen::VectorXd f;
    if(std::all_of(bounds_cond.cbegin(), bounds_cond.cend(), [](const std::tuple<boundary_type, std::function<double(double, double)>> &bound)
                                                             { return std::get<boundary_type>(bound) == boundary_type::FLOW; }))
    {
        K.resize(mesh.nodes_count()+1, mesh.nodes_count()+1);
        f = Eigen::VectorXd::Zero(mesh.nodes_count()+1);
        f[mesh.nodes_count()] = volume;
    }
    else 
    {
        K.resize(mesh.nodes_count(), mesh.nodes_count());
        K_bound.resize(mesh.nodes_count(), mesh.nodes_count());
        f = Eigen::VectorXd::Zero(mesh.nodes_count());
    }

    double time = omp_get_wtime();
    integrate_right_part(mesh, f, right_part);
    std::cout << "Right part Integrate: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    create_matrix(mesh, bounds_cond, K, K_bound, integrate_gradient_pair<double>, p1, influence_fun);
    if(std::all_of(bounds_cond.cbegin(), bounds_cond.cend(), [](const std::tuple<boundary_type, std::function<double(double, double)>> &bound)
                                                             { return std::get<boundary_type>(bound) == boundary_type::FLOW; }))
        K += nonlocal_condition(mesh);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;
    std::cout << "Nonzero elements count = " << K.nonZeros() + K_bound.nonZeros() << std::endl;

    time = omp_get_wtime();
    boundary_condition(mesh, temperature_nodes_vectors(mesh, bounds_cond), bounds_cond, 1., K_bound, f);
    std::cout << "Boundary filling: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::VectorXd x = solver.solve(f);
    std::cout << "System solving: " << omp_get_wtime() - time << std::endl;
    
    //save_as_vtk(path, mesh, x);
    mesh.print_to_file(path, x.data());
}

void nonstationary(const std::string &path,
                   const mesh_2d<double> &mesh, const double tau, const size_t time_steps,
                   const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                   const std::function<double(double, double)> &init_dist,
                   const std::function<double(double, double)> &right_part,
                   const double p1, const std::function<double(double, double, double, double)> &influence_fun,
                   const uint64_t print_frequency)
{
    Eigen::SparseMatrix<double> K      (mesh.nodes_count(), mesh.nodes_count()),
                                K_bound(mesh.nodes_count(), mesh.nodes_count()),
                                C      (mesh.nodes_count(), mesh.nodes_count()),
                                C_bound(mesh.nodes_count(), mesh.nodes_count());
    create_matrix(mesh, bounds_cond, K, K_bound, integrate_gradient_pair<double>, p1, influence_fun);
    create_matrix(mesh, bounds_cond, C, C_bound, integrate_basic_pair<double>, 1., influence_fun);
    C_bound.setZero();
    K *= tau;
    K_bound *= tau;
    K += C;

    std::vector<std::vector<uint32_t>> temperature_nodes = temperature_nodes_vectors(mesh, bounds_cond);
    for(size_t b = 0; b < temperature_nodes.size(); ++b)
        for(auto node : temperature_nodes[b])
            K.coeffRef(node, node) = 1.;

    Eigen::VectorXd f = Eigen::VectorXd::Zero(mesh.nodes_count()),
                    T_prev(mesh.nodes_count()), T(mesh.nodes_count());
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        T_prev[i] = init_dist(mesh.coord(i, 0), mesh.coord(i, 1));

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);

    if(print_frequency != uint64_t(-1))
    {
        mesh.print_to_file(path+std::string("0.csv"), T_prev.data());
        std::cout << "step = " << 0 << " Volume = " << integrate_solution(mesh, T_prev) << std::endl;
    }
    for(size_t i = 1; i < time_steps; ++i)
    {
        f = C.selfadjointView<Eigen::Lower>() * T_prev;
        integrate_right_part(mesh, f, right_part);
        boundary_condition(mesh, temperature_nodes, bounds_cond, tau, K_bound, f);
        T = solver.solve(f);
        if(i % print_frequency == 0)
        {
            mesh.print_to_file(path+std::to_string(i)+std::string(".csv"), T.data());
            std::cout << "step = " << i << " Volume = " << integrate_solution(mesh, T) << std::endl;
        }
        T_prev.swap(T);
    }
}

}