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

#pragma GCC diagnostic pop
};

static std::tuple<std::vector<Eigen::Triplet<double, node_info>>, size_t, std::vector<Eigen::Triplet<double, node_info>>, size_t>
    mesh_analysis(const mesh_2d<double> &mesh, const std::set<node_info> &fixed_nodes)
{
    std::vector<uint32_t> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0);

    std::function<void(node_info, node_info, size_t)>
        counter_loc = [&mesh, &fixed_nodes, &shifts_loc, &shifts_bound_loc]
                      (node_info node_i, node_info node_j, size_t el)
                      {
                          if(fixed_nodes.find({mesh.node(el, node_i), node_i.comp}) == fixed_nodes.cend() &&
                             fixed_nodes.find({mesh.node(el, node_j), node_j.comp}) == fixed_nodes.cend())
                              ++shifts_loc[el+1];
                          else if(mesh.node(el, node_i) != mesh.node(el, node_j))
                              ++shifts_bound_loc[el+1];
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
                         if(fixed_nodes.find({mesh.node(el, node_i), node_i.comp}) == fixed_nodes.cend() &&
                            fixed_nodes.find({mesh.node(el, node_j), node_j.comp}) == fixed_nodes.cend())
                             triplets[shifts_loc[el]++] = Eigen::Triplet<double, node_info>(node_i, node_j, *reinterpret_cast<double*>(&el));
                         else if(mesh.node(el, node_i) != mesh.node(el, node_j))
                             triplets_bound[shifts_bound_loc[el]++] = fixed_nodes.find({mesh.node(el, node_j), node_j.comp}) == fixed_nodes.cend() ?
                                                                      Eigen::Triplet<double, node_info>(node_j, node_i, *reinterpret_cast<double*>(&el)) :
                                                                      Eigen::Triplet<double, node_info>(node_i, node_j, *reinterpret_cast<double*>(&el));
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

template<class Type>
static Type integrate(const finite_element::element_2d_integrate_base<Type> *e,
                      const node_info i, const node_info j, const matrix<Type> &jacobi_matrices)
{
    static constexpr Type nu = 0.3, E = 2.1e5;
    static constexpr Type A = E / (1-nu*nu), B = nu * A, C = 0.5 * E / (1. + nu);

    Type coeff_1 = i.comp == j.comp ? A : B, coeff_2 = C;
    Type integral = 0.0;
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        integral += (coeff_1 * (e->qNxi(i, q)*jacobi_matrices(q, 3) + e->qNeta(i, q)*jacobi_matrices(q, 2)) *
                               (e->qNxi(j, q)*jacobi_matrices(q, 3) + e->qNeta(j, q)*jacobi_matrices(q, 2)) +
                     coeff_2 * (e->qNxi(i, q)*jacobi_matrices(q, 1) + e->qNeta(i, q)*jacobi_matrices(q, 0)) *
                               (e->qNxi(j, q)*jacobi_matrices(q, 1) + e->qNeta(j, q)*jacobi_matrices(q, 0))) /
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2)) * e->weight(q);
    return integral;
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
        triplets[i] = Eigen::Triplet<double, node_info>({2*mesh.node(el, triplets[i].row())+triplets[i].row().comp, triplets[i].row().comp},
                                                        {2*mesh.node(el, triplets[i].col())+triplets[i].col().comp, triplets[i].col().comp},
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

void create_matrix(const mesh_2d<double> &mesh,
                          const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                       boundary_type, std::function<double(double, double)>>> &bounds_cond/*,
                          Eigen::SparseMatrix<double> &K, Eigen::SparseMatrix<double> &K_bound*/)
{
    std::set<node_info> fixed_nodes = fixed_nodes_set(mesh, bounds_cond);
    auto [triplets, classic_count, triplets_bound, classic_bound_count] = mesh_analysis(mesh, fixed_nodes);

    
    triplets_run_loc(mesh, triplets,       fixed_nodes.size(), classic_count      );
    triplets_run_loc(mesh, triplets_bound,                  0, classic_bound_count);

/*
    for(auto it = fixed_nodes.cbegin(); it != fixed_nodes.cend(); ++it)
        std::cout << it->comp << " " << it->number << std::endl;

    std::cout << "=======================" << std::endl;

    for(auto it = triplets.cbegin(); it != triplets.cend(); ++it)
        std::cout << (it->col().comp ? "Y" : "X") << " " << it->col().number << " "
                  << (it->row().comp ? "Y" : "X") << " " << it->row().number << " "
                  << *reinterpret_cast<const double*>(&it->value()) << std::endl;

    std::cout << "=======================" << std::endl;

    for(auto it = triplets_bound.cbegin(); it != triplets_bound.cend(); ++it)
        std::cout << (it->col().comp ? "Y" : "X") << " " << it->col().number << " "
                  << (it->row().comp ? "Y" : "X") << " " << it->row().number << " "
                  << *reinterpret_cast<const double*>(&it->value()) << std::endl;

    std::cout << triplets.size() << " " << triplets_bound.size() << std::endl;
*/

    Eigen::SparseMatrix<double> K(2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    
    fixed_nodes.clear();
    K_bound.setFromTriplets(triplets_bound.begin(), triplets_bound.end());
    triplets_bound.clear();
    K.setFromTriplets(triplets.cbegin(), triplets.cend());

    std::cout << Eigen::MatrixXd(K) << std::endl << std::endl
              << Eigen::MatrixXd(K_bound) << std::endl;
}

/*
void stationary(const std::string &path, const mesh_2d<double> &mesh,
                )
                */