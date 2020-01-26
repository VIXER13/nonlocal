#include <set>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstring>

#include "omp.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "finite_element_routine.hpp"

#include "static_analysis.hpp"

struct node_info {
    uint32_t              number : 31;
    enum component {X, Y} comp   :  1;

    node_info(uint64_t number, component comp = X) :
        number(static_cast<uint32_t>(number)), comp(comp) {}
    
    operator Eigen::SparseMatrix<double>::StorageIndex() const {
        return static_cast<Eigen::SparseMatrix<double>::StorageIndex>(number);
    }

    template<class Type>
    friend bool operator<(const node_info left, const Type right) {
        return left.number < right;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"

    friend bool operator<(const node_info left, const node_info right) {
        return *reinterpret_cast<const uint32_t*>(&left) < *reinterpret_cast<const uint32_t*>(&right);
    }

    friend bool operator!=(const node_info left, const node_info right) {
        return *reinterpret_cast<const uint32_t*>(&left) != *reinterpret_cast<const uint32_t*>(&right);
    }

#pragma GCC diagnostic pop
};

static void save_as_vtk(const std::string &path, const mesh_2d<double> &mesh, const Eigen::VectorXd &u)
{
    std::ofstream fout(path);
    fout.precision(20);

    fout << "# vtk DataFile Version 4.2" << std::endl
         << "Displacement"               << std::endl
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

    fout << "SCALARS U_X double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < static_cast<size_t>(u.size() / 2); ++i)
        fout << u[2*i] << std::endl;

    fout << "SCALARS U_Y double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < static_cast<size_t>(u.size() / 2); ++i)
        fout << u[2*i+1] << std::endl;
}

template<class Type>
static Type integrate(const finite_element::element_2d_integrate_base<Type> *e,
                      const node_info i, const node_info j, const matrix<Type> &jacobi_matrices)
{
    static constexpr Type nu = 0.3, E = 2.1e5;
    static constexpr Type A = E / (1-nu*nu), B = nu * A, C = 0.5 * E / (1. + nu);

    Type coeff_1 = i.comp == j.comp ? A : B, coeff_2 = C;
    Type integral = 0.0;
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        integral += (coeff_1 * ( e->qNxi(i, q)*jacobi_matrices(q, 3) - e->qNeta(i, q)*jacobi_matrices(q, 2)) *
                               ( e->qNxi(j, q)*jacobi_matrices(q, 3) - e->qNeta(j, q)*jacobi_matrices(q, 2)) +
                     coeff_2 * (-e->qNxi(i, q)*jacobi_matrices(q, 1) + e->qNeta(i, q)*jacobi_matrices(q, 0)) *
                               (-e->qNxi(j, q)*jacobi_matrices(q, 1) + e->qNeta(j, q)*jacobi_matrices(q, 0))) /
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2)) * e->weight(q);
    return integral;
}

template<class Type>
static Type integrate_force_bound(const finite_element::element_1d_integrate_base<Type> *be, const size_t i,
                                  const matrix<Type> &coords, const matrix<Type> &jacobi_matrices, 
                                  const std::function<Type(Type, Type)> &fun)
{
    Type integral = 0.0;
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        integral += fun(coords(q, 0), coords(q, 1)) * be->weight(q) * be->qN(i, q) *
                    sqrt(jacobi_matrices(q, 0)*jacobi_matrices(q, 0) + jacobi_matrices(q, 1)*jacobi_matrices(q, 1));
    return integral;
}

static std::tuple<std::vector<Eigen::Triplet<double, node_info>>, size_t, std::vector<Eigen::Triplet<double, node_info>>, size_t>
    mesh_analysis(const mesh_2d<double> &mesh, const std::set<node_info> &fixed_nodes)
{
    std::vector<uint32_t> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0);

    std::function<void(node_info, node_info, size_t)>
        counter_loc = [&mesh, &fixed_nodes, &shifts_loc, &shifts_bound_loc]
                      (node_info node_i, node_info node_j, size_t el)
                      {
                          node_info glob_i = {mesh.node_number(el, node_i.number), node_i.comp},
                                    glob_j = {mesh.node_number(el, node_j.number), node_j.comp};
                          if(glob_i >= glob_j)
                          {
                            if(fixed_nodes.find(glob_i) == fixed_nodes.cend() &&
                                fixed_nodes.find(glob_j) == fixed_nodes.cend())
                                ++shifts_loc[el+1];
                            else if(glob_i != glob_j)
                                ++shifts_bound_loc[el+1];
                          }
                      };
    mesh_run_loc(mesh, [&counter_loc](size_t i, size_t j, size_t el)
                       {
                           counter_loc({i, node_info::X}, {j, node_info::X}, el);
                           counter_loc({i, node_info::X}, {j, node_info::Y}, el);
                           counter_loc({i, node_info::Y}, {j, node_info::X}, el);
                           counter_loc({i, node_info::Y}, {j, node_info::Y}, el);
                       }
                );

    shifts_loc[0] = fixed_nodes.size();
    for(size_t i = 1; i < shifts_loc.size(); ++i)
    {
        shifts_loc[i] += shifts_loc[i-1];
        shifts_bound_loc[i] += shifts_bound_loc[i-1];
    }

    std::vector<Eigen::Triplet<double, node_info>> triplets(shifts_loc.back()),
                                                   triplets_bound(shifts_bound_loc.back());
    for(auto [it, i] = std::make_tuple(fixed_nodes.cbegin(), static_cast<size_t>(0)); it != fixed_nodes.cend(); ++it, ++i)
    {
        uint32_t index = 2 * it->number + static_cast<uint32_t>(it->comp);
        triplets[i] = Eigen::Triplet<double, node_info>({index, it->comp}, {index, it->comp}, 1.);
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"

    std::function<void(node_info, node_info, size_t)>
        setter_loc = [&mesh, &fixed_nodes, &shifts_loc, &shifts_bound_loc, &triplets, &triplets_bound]
                     (node_info node_i, node_info node_j, size_t el)
                     {
                         node_info glob_i = {mesh.node_number(el, node_i.number), node_i.comp},
                                   glob_j = {mesh.node_number(el, node_j.number), node_j.comp};
                        if(glob_i >= glob_j)
                        {
                            if(fixed_nodes.find(glob_i) == fixed_nodes.cend() &&
                               fixed_nodes.find(glob_j) == fixed_nodes.cend())
                                triplets[shifts_loc[el]++] = Eigen::Triplet<double, node_info>(node_i, node_j, *reinterpret_cast<double*>(&el));
                            else if(glob_i != glob_j)
                                triplets_bound[shifts_bound_loc[el]++] = fixed_nodes.find(glob_j) == fixed_nodes.cend() ?
                                                                        Eigen::Triplet<double, node_info>(node_j, node_i, *reinterpret_cast<double*>(&el)) :
                                                                        Eigen::Triplet<double, node_info>(node_i, node_j, *reinterpret_cast<double*>(&el));
                        }
                     };

#pragma GCC diagnostic pop

    mesh_run_loc(mesh, [&setter_loc](size_t i, size_t j, size_t el)
                       {
                           setter_loc({i, node_info::X}, {j, node_info::X}, el);
                           setter_loc({i, node_info::X}, {j, node_info::Y}, el);
                           setter_loc({i, node_info::Y}, {j, node_info::X}, el);
                           setter_loc({i, node_info::Y}, {j, node_info::Y}, el);
                       }
                );
    
    return {std::move(triplets), shifts_loc.back(), std::move(triplets_bound), shifts_bound_loc.back()};
}

static void triplets_run_loc(const mesh_2d<double> &mesh, std::vector<Eigen::Triplet<double, node_info>> &triplets,
                             const size_t start, const size_t finish)
{
#pragma omp parallel default(none) shared(mesh, triplets)
{
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    size_t el = 0, el_prev = size_t(-1);
    matrix<double> jacobi_matrices;
#pragma omp for
    for(size_t i = start; i < finish; ++i)
    {
        el = *reinterpret_cast<const size_t*>(&triplets[i].value());
        if(el != el_prev && mesh.element_type(el) != mesh.element_type(el_prev))
        {
            e = mesh.element_2d(mesh.element_type(el));
            approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        }
        el_prev = el;
        triplets[i] = Eigen::Triplet<double, node_info>({2*mesh.node_number(el, triplets[i].row())+triplets[i].row().comp, triplets[i].row().comp},
                                                        {2*mesh.node_number(el, triplets[i].col())+triplets[i].col().comp, triplets[i].col().comp},
                                                        integrate(e, triplets[i].row(), triplets[i].col(), jacobi_matrices));
    }
}
}

static std::set<node_info> fixed_nodes_set(const mesh_2d<double> &mesh,
                                           const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                                        boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::set<node_info> fixed_nodes;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                fixed_nodes.insert(node_info(*node, node_info::X));

        if(std::get<2>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                fixed_nodes.insert(node_info(*node, node_info::Y));
    }
    return fixed_nodes;
}

static void create_matrix(const mesh_2d<double> &mesh,
                          const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                       boundary_type, std::function<double(double, double)>>> &bounds_cond,
                          Eigen::SparseMatrix<double> &K, Eigen::SparseMatrix<double> &K_bound)
{
    std::set<node_info> fixed_nodes = fixed_nodes_set(mesh, bounds_cond);
    auto [triplets, classic_count, triplets_bound, classic_bound_count] = mesh_analysis(mesh, fixed_nodes);

    triplets_run_loc(mesh, triplets,       fixed_nodes.size(), classic_count      );
    triplets_run_loc(mesh, triplets_bound,                  0, classic_bound_count);

    fixed_nodes.clear();
    K_bound.setFromTriplets(triplets_bound.begin(), triplets_bound.end());
    triplets_bound.clear();
    K.setFromTriplets(triplets.cbegin(), triplets.cend());

    //std::cout << Eigen::MatrixXd(K) << std::endl << std::endl
    //          << Eigen::MatrixXd(K_bound) << std::endl;
}

static void translation(const mesh_2d<double> &mesh, const Eigen::SparseMatrix<double> &K_bound, Eigen::VectorXd &f,
                        const std::function<double(double, double)> &boundaryFun, const size_t node)
{
    double temp = boundaryFun(mesh.coord(node >> 1, 0), mesh.coord(node >> 1, 1));
    for(typename Eigen::SparseMatrix<double>::InnerIterator it(K_bound, node); it; ++it)
        f[it.row()] -= temp * it.value();
}

static void boundary_condition(const mesh_2d<double> &mesh, const std::vector<std::vector<uint32_t>> &temperature_nodes,
                               const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                            boundary_type, std::function<double(double, double)>>> &bounds_cond,
                               const double tau, const Eigen::SparseMatrix<double> &K_bound, Eigen::VectorXd &f)
{
    const finite_element::element_1d_integrate_base<double> *be = nullptr;
    matrix<double> coords, jacobi_matrices;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::FORCE)
            for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
            {
                be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
                approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
                approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
                for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
                    f[mesh.boundary(b)(el, i)] += tau*integrate_force_bound(be, i, coords, jacobi_matrices, std::get<1>(bounds_cond[b]));
            }

        if(std::get<2>(bounds_cond[b]) == boundary_type::FORCE)
            for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
            {
                be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
                approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
                approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
                for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
                    f[mesh.boundary(b)(el, i)] += tau*integrate_force_bound(be, i, coords, jacobi_matrices, std::get<3>(bounds_cond[b]));
            }
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node : temperature_nodes[b])
                translation(mesh, K_bound, f, std::get<1>(bounds_cond[b]), 2*node);

        if(std::get<2>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node : temperature_nodes[b])
                translation(mesh, K_bound, f, std::get<3>(bounds_cond[b]), 2*node+1);
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node : temperature_nodes[b])
                f[2*node]   = std::get<1>(bounds_cond[b])(mesh.coord(node, 0), mesh.coord(node, 1));

        if(std::get<2>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto node : temperature_nodes[b])
                f[2*node+1] = std::get<3>(bounds_cond[b])(mesh.coord(node, 0), mesh.coord(node, 1));
    }
}

static std::vector<std::vector<uint32_t>> fixed_nodes_vectors(const mesh_2d<double> &mesh,
                                              const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                                           boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::vector<std::vector<uint32_t>> fixed_nodes(bounds_cond.size());
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(std::get<0>(bounds_cond[b]) == boundary_type::FIXED ||
           std::get<2>(bounds_cond[b]) == boundary_type::FIXED)
            for(auto [node, k] = std::make_tuple(mesh.boundary(b).cbegin(), static_cast<size_t>(0)); node != mesh.boundary(b).cend(); ++node)
            {
                for(k = 0; k < fixed_nodes.size(); ++k)
                    if(std::find(fixed_nodes[k].cbegin(), fixed_nodes[k].cend(), *node) != fixed_nodes[k].cend())
                        break;
                if(k == fixed_nodes.size())
                    fixed_nodes[b].push_back(*node);
            }
    return fixed_nodes;
}

void strains_calculate(const mesh_2d<double> &mesh, const Eigen::VectorXd &u,
                       Eigen::VectorXd &eps11, Eigen::VectorXd &eps12, Eigen::VectorXd &eps22)
{
    

    std::vector<uint8_t> count(mesh.nodes_count());
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ++count[mesh.node_number(el, i)];
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        eps11[i] /= count[i];
        eps12[i] /= count[i];
        eps22[i] /= count[i];
    }
}

void stationary(const std::string &path, const mesh_2d<double> &mesh,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                             boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    Eigen::VectorXd f = Eigen::VectorXd::Zero(2*mesh.nodes_count());;
    Eigen::SparseMatrix<double> K(2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    
    create_matrix(mesh, bounds_cond, K, K_bound);

    boundary_condition(mesh, fixed_nodes_vectors(mesh, bounds_cond), bounds_cond, 1., K_bound, f);

    //std::cout << f.transpose() << std::endl << std::endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::VectorXd u = solver.solve(f);

    Eigen::VectorXd eps11 = Eigen::VectorXd::Zero(mesh.nodes_count()),
                    eps12 = Eigen::VectorXd::Zero(mesh.nodes_count()),
                    eps22 = Eigen::VectorXd::Zero(mesh.nodes_count());

    save_as_vtk(path, mesh, u);

    //std::cout << x.transpose() << std::endl;

    std::ofstream fout_x(std::string("results//text_x.csv")),
                  fout_y(std::string("results//text_y.csv"));
    fout_x.precision(20);
    fout_y.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        fout_x << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i) << std::endl;
        fout_y << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i+1) << std::endl;
    }
}