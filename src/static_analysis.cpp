#include <set>
#include <algorithm>
#include "omp.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/PardisoSupport"
#include "finite_element_routine.hpp"
#include "static_analysis.hpp"

namespace statics_with_nonloc
{

// Небольшая структура, которая объединяет в себе номер узла и индекс переменной, который на ней задан.
struct node_info
{
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

static void save_as_vtk(const std::string &path,            const mesh_2d<double> &mesh,        const Eigen::VectorXd &u,
                        const std::vector<double> &eps11,   const std::vector<double> &eps22,   const std::vector<double> &eps12,
                        const std::vector<double> &sigma11, const std::vector<double> &sigma22, const std::vector<double> &sigma12)
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

    fout << "SCALARS U_X double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < size_t(u.size() / 2); ++i)
        fout << u[2*i] << std::endl;

    fout << "SCALARS U_Y double " << 1 << std::endl
         << "LOOKUP_TABLE default" << std::endl;
    for(size_t i = 0; i < size_t(u.size() / 2); ++i)
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

template<class Type>
static Type integrate_loc(const finite_element::element_2d_integrate_base<Type> *const e,
                          const node_info i, const node_info j, const matrix<Type> &jacobi_matrices, size_t shift,
                          const std::array<Type, 3> &coeff)
{
    Type integral = 0.;
    if(i.comp == node_info::X)
    {
        if(j.comp == node_info::X) // XX
        {
            for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
                integral += (coeff[0] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                        ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2)) +
                             coeff[2] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                        (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0))) /
                            (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
        }
        else // XY
        {
            for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
                integral += (coeff[1] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                        ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2)) +
                             coeff[2] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                        (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0))) /
                            (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
        }
    }
    else
    {
        if(j.comp == node_info::X) //YX
        {
            for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
                integral += (coeff[1] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                        (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0)) +
                             coeff[2] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                        ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2))) /
                            (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
        }
        else // YY
        {
            for(size_t q = 0; q < e->qnodes_count(); ++q, ++shift)
                integral += (coeff[0] * (-e->qNxi(i, q)*jacobi_matrices(shift, 1) + e->qNeta(i, q)*jacobi_matrices(shift, 0)) *
                                        (-e->qNxi(j, q)*jacobi_matrices(shift, 1) + e->qNeta(j, q)*jacobi_matrices(shift, 0)) +
                             coeff[2] * ( e->qNxi(i, q)*jacobi_matrices(shift, 3) - e->qNeta(i, q)*jacobi_matrices(shift, 2)) *
                                        ( e->qNxi(j, q)*jacobi_matrices(shift, 3) - e->qNeta(j, q)*jacobi_matrices(shift, 2))) /
                            (jacobi_matrices(shift, 0)*jacobi_matrices(shift, 3) - jacobi_matrices(shift, 1)*jacobi_matrices(shift, 2)) * e->weight(q);
        }
    }   
    return integral;
}

template<class Type>
static Type integrate_nonloc(const finite_element::element_2d_integrate_base<Type> *const eL,
                             const finite_element::element_2d_integrate_base<Type> *const eNL,
                             const node_info iL, const node_info jNL, size_t shiftL, size_t shiftNL,
                             const matrix<Type> &coords, const matrix<Type> &jacobi_matrices,
                             const std::function<Type(Type, Type, Type, Type)> &influence_fun,
                             const std::array<Type, 3> &coeff)
{
    const size_t sub = shiftNL;
    Type integral = 0.;
    if(iL.comp == node_info::X)
    {
        if(jNL.comp == node_info::X) // XX
        {
            Type int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
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
                            (coeff[0] * int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)) +
                             coeff[2] * int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)));                
            }
        }
        else // XY
        {
            Type int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
            for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL)
            {
                int_with_weight_x = 0.;
                int_with_weight_y = 0.;
                for(size_t qNL = 0, shiftNL = sub; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL)
                {
                    finit = eNL->weight(qNL) * influence_fun(coords(shiftL, 0), coords(shiftNL, 0), coords(shiftL, 1), coords(shiftNL, 1));
                    int_with_weight_x += finit * ( eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 1) - eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 0));
                    int_with_weight_y += finit * (-eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 3) + eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 2));
                }
                integral += eL->weight(qL) *
                            (coeff[1] * int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)) +
                             coeff[2] * int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)));                
            }
        }
    }
    else
    {
        if(jNL.comp == node_info::X) //YX
        {
            Type int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
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
                            (coeff[1] * int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)) +
                             coeff[2] * int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)));                
            }
        }
        else // YY
        {
            Type int_with_weight_x = 0., int_with_weight_y = 0., finit = 0.;
            for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL)
            {
                int_with_weight_x = 0.;
                int_with_weight_y = 0.;
                for(size_t qNL = 0, shiftNL = sub; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL)
                {
                    finit = eNL->weight(qNL) * influence_fun(coords(shiftL, 0), coords(shiftNL, 0), coords(shiftL, 1), coords(shiftNL, 1));
                    int_with_weight_x += finit * ( eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 1) - eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 0));
                    int_with_weight_y += finit * (-eNL->qNxi(jNL, qNL) * jacobi_matrices(shiftNL, 3) + eNL->qNeta(jNL, qNL) * jacobi_matrices(shiftNL, 2));
                }
                integral += eL->weight(qL) *
                            (coeff[0] * int_with_weight_x * ( eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 1) - eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 0)) +
                             coeff[2] * int_with_weight_y * (-eL->qNxi(iL, qL) * jacobi_matrices(shiftL, 3) + eL->qNeta(iL, qL) * jacobi_matrices(shiftL, 2)));                
            }
        }
    }   
    return integral;
}

template<class Type>
static Type integrate_force_bound(const finite_element::element_1d_integrate_base<Type> *be, const size_t i,
                                  const matrix<Type> &coords, const matrix<Type> &jacobi_matrices, 
                                  const std::function<Type(Type, Type)> &fun)
{
    Type integral = 0.;
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        integral += fun(coords(q, 0), coords(q, 1)) * be->weight(q) * be->qN(i, q) *
                    sqrt(jacobi_matrices(q, 0)*jacobi_matrices(q, 0) + jacobi_matrices(q, 1)*jacobi_matrices(q, 1));
    return integral;
}

static std::set<node_info> kinematic_nodes_set(const mesh_2d<double> &mesh,
                                               const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                                            boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::set<node_info> kinematic_nodes;
    for(size_t b = 0; b < bounds_cond.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                kinematic_nodes.insert(node_info(*node, node_info::X));

        if(std::get<2>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node = mesh.boundary(b).cbegin(); node != mesh.boundary(b).cend(); ++node)
                kinematic_nodes.insert(node_info(*node, node_info::Y));
    }
    return std::move(kinematic_nodes);
}

static std::vector<std::vector<uint32_t>> kinematic_nodes_vectors(const mesh_2d<double> &mesh,
                                              const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                                           boundary_type, std::function<double(double, double)>>> &bounds_cond)
{
    std::vector<std::vector<uint32_t>> kinematic_nodes(bounds_cond.size());
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        if(std::get<0>(bounds_cond[b]) == boundary_type::TRANSLATION ||
           std::get<2>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto [node, k] = std::make_tuple(mesh.boundary(b).cbegin(), static_cast<size_t>(0)); node != mesh.boundary(b).cend(); ++node)
            {
                for(k = 0; k < kinematic_nodes.size(); ++k)
                    if(std::find(kinematic_nodes[k].cbegin(), kinematic_nodes[k].cend(), *node) != kinematic_nodes[k].cend())
                        break;
                if(k == kinematic_nodes.size())
                    kinematic_nodes[b].push_back(*node);
            }
    return std::move(kinematic_nodes);
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
                    f[2*mesh.boundary(b)(el, i)] += tau*integrate_force_bound(be, i, coords, jacobi_matrices, std::get<1>(bounds_cond[b]));
            }

        if(std::get<2>(bounds_cond[b]) == boundary_type::FORCE)
            for(size_t el = 0; el < mesh.boundary(b).rows(); ++el)
            {
                be = mesh.element_1d(mesh.elements_on_bound_types(b)[el]);
                approx_jacobi_matrices_bound(mesh, be, b, el, jacobi_matrices);
                approx_quad_nodes_coord_bound(mesh, be, b, el, coords);
                for(size_t i = 0; i < mesh.boundary(b).cols(el); ++i)
                    f[2*mesh.boundary(b)(el, i)+1] += tau*integrate_force_bound(be, i, coords, jacobi_matrices, std::get<3>(bounds_cond[b]));
            }
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node : temperature_nodes[b])
                translation(mesh, K_bound, f, std::get<1>(bounds_cond[b]), 2*node);

        if(std::get<2>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node : temperature_nodes[b])
                translation(mesh, K_bound, f, std::get<3>(bounds_cond[b]), 2*node+1);
    }

    for(size_t b = 0; b < temperature_nodes.size(); ++b)
    {
        if(std::get<0>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node : temperature_nodes[b])
                f[2*node]   = std::get<1>(bounds_cond[b])(mesh.coord(node, 0), mesh.coord(node, 1));

        if(std::get<2>(bounds_cond[b]) == boundary_type::TRANSLATION)
            for(auto node : temperature_nodes[b])
                f[2*node+1] = std::get<3>(bounds_cond[b])(mesh.coord(node, 0), mesh.coord(node, 1));
    }
}

static std::array<std::vector<uint32_t>, 4>
    mesh_analysis(const mesh_2d<double> &mesh, const std::set<node_info> &kinematic_nodes, const bool nonlocal)
{
    std::vector<uint32_t> shifts_loc(mesh.elements_count()+1, 0), shifts_bound_loc(mesh.elements_count()+1, 0),
                          shifts_nonloc, shifts_bound_nonloc;

    const auto counter_loc = 
        [&mesh, &kinematic_nodes, &shifts_loc, &shifts_bound_loc]
        (node_info node_i, node_info node_j, size_t el)
        {
            node_info glob_i = {mesh.node_number(el, node_i.number), node_i.comp},
                      glob_j = {mesh.node_number(el, node_j.number), node_j.comp};
            if(2 * glob_i.number + glob_i.comp >= 2 * glob_j.number + glob_j.comp)
            {
                if(kinematic_nodes.find(glob_i) == kinematic_nodes.cend() &&
                   kinematic_nodes.find(glob_j) == kinematic_nodes.cend())
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
                       });

    shifts_loc[0] = kinematic_nodes.size();
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
            [&mesh, &kinematic_nodes, &shifts_nonloc, &shifts_bound_nonloc]
            (node_info node_iL, node_info node_jNL, size_t elL, size_t elNL)
            {
                node_info glob_iL  = {mesh.node_number(elL,  node_iL.number),  node_iL.comp},
                          glob_jNL = {mesh.node_number(elNL, node_jNL.number), node_jNL.comp};
                if(2 * glob_iL.number + glob_iL.comp >= 2 * glob_jNL.number + glob_jNL.comp)
                {
                    if(kinematic_nodes.find(glob_iL)  == kinematic_nodes.cend() &&
                       kinematic_nodes.find(glob_jNL) == kinematic_nodes.cend())
                        ++shifts_nonloc[elL+1];
                    else if(glob_iL != glob_jNL)
                        ++shifts_bound_nonloc[elL+1];
                }
            };

        mesh_run_nonloc(mesh, [&counter_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
                              {
                                  counter_nonloc({iL, node_info::X}, {jNL, node_info::X}, elL, elNL);
                                  counter_nonloc({iL, node_info::X}, {jNL, node_info::Y}, elL, elNL);
                                  counter_nonloc({iL, node_info::Y}, {jNL, node_info::X}, elL, elNL);
                                  counter_nonloc({iL, node_info::Y}, {jNL, node_info::Y}, elL, elNL);
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

static std::array<std::vector<Eigen::Triplet<double, node_info>>, 2>
    triplets_fill(const mesh_2d<double> &mesh, const std::set<node_info> &kinematic_nodes, const parameters<double> &params,
                  const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    static constexpr double MAX_LOCAL_WEIGHT = 0.999;
    bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
    auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, kinematic_nodes, nonlocal);
    std::vector<Eigen::Triplet<double, node_info>> triplets      (nonlocal ? shifts_nonloc.back()       : shifts_loc.back()),
                                                   triplets_bound(nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
    std::cout << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;
    for(auto [it, i] = std::make_tuple(kinematic_nodes.cbegin(), size_t(0)); it != kinematic_nodes.cend(); ++it, ++i)
    {
        uint32_t index = 2 * it->number + it->comp;
        triplets[i] = Eigen::Triplet<double, node_info>({index, it->comp}, {index, it->comp}, 1.);
    }

    const std::vector<uint32_t> shifts_quad = quadrature_shifts_init(mesh);
    const matrix<double> all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts_quad);
    const std::array<double, 3> coeffs = {            params.E / (1. - params.nu*params.nu),
                                          params.nu * params.E / (1. - params.nu*params.nu),
                                                0.5 * params.E / (1. + params.nu)           };

    const auto filler_loc =
        [&mesh, &kinematic_nodes, &shifts_loc, &shifts_bound_loc, &triplets, &triplets_bound, &shifts_quad, &all_jacobi_matrices, p1, &coeffs]
        (node_info node_i, node_info node_j, size_t el)
        {
            node_info glob_i = {mesh.node_number(el, node_i.number), node_i.comp},
                      glob_j = {mesh.node_number(el, node_j.number), node_j.comp},
                      row = 2 * glob_i.number + glob_i.comp,
                      col = 2 * glob_j.number + glob_j.comp;
            if(row >= col)
            {
                double integral = p1 * integrate_loc(mesh.element_2d(mesh.element_type(el)), node_i, node_j, all_jacobi_matrices, shifts_quad[el], coeffs);
                if(kinematic_nodes.find(glob_i) == kinematic_nodes.cend() &&
                   kinematic_nodes.find(glob_j) == kinematic_nodes.cend())
                    triplets[shifts_loc[el]++] = Eigen::Triplet<double, node_info>(row, col, integral);
                else if(glob_i != glob_j)
                    triplets_bound[shifts_bound_loc[el]++] = kinematic_nodes.find(glob_j) == kinematic_nodes.cend() ?
                                                             Eigen::Triplet<double, node_info>(col, row, integral) :
                                                             Eigen::Triplet<double, node_info>(row, col, integral);
            }
        };

    mesh_run_loc(mesh, [&filler_loc](size_t i, size_t j, size_t el)
                       {
                           filler_loc({i, node_info::X}, {j, node_info::X}, el);
                           filler_loc({i, node_info::X}, {j, node_info::Y}, el);
                           filler_loc({i, node_info::Y}, {j, node_info::X}, el);
                           filler_loc({i, node_info::Y}, {j, node_info::Y}, el);
                       });

    if(nonlocal)
    {
        const matrix<double> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts_quad);

        const auto filler_nonloc =
            [&mesh, &kinematic_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
             &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1, &coeffs]
            (node_info node_iL, node_info node_jNL, size_t elL, size_t elNL)
            {
                node_info glob_iL  = {mesh.node_number(elL,  node_iL.number),  node_iL.comp},
                          glob_jNL = {mesh.node_number(elNL, node_jNL.number), node_jNL.comp},
                          row = 2 * glob_iL.number  + glob_iL.comp,
                          col = 2 * glob_jNL.number + glob_jNL.comp;
                if(row >= col)
                {
                    double integral = p2 * integrate_nonloc(mesh.element_2d(mesh.element_type(elL )),
                                                            mesh.element_2d(mesh.element_type(elNL)), 
                                                            node_iL, node_jNL, shifts_quad[elL], shifts_quad[elNL],
                                                            all_quad_coords, all_jacobi_matrices, influence_fun, coeffs);
                    if(kinematic_nodes.find(glob_iL)  == kinematic_nodes.cend() &&
                       kinematic_nodes.find(glob_jNL) == kinematic_nodes.cend())
                        triplets[shifts_nonloc[elL]++] = Eigen::Triplet<double, node_info>(row, col, integral);
                    else if(glob_iL != glob_jNL)
                        triplets_bound[shifts_bound_nonloc[elL]++] = kinematic_nodes.find(glob_jNL) == kinematic_nodes.cend() ?
                                                                     Eigen::Triplet<double, node_info>(col, row,  integral) :
                                                                     Eigen::Triplet<double, node_info>(row, col, integral);
                }
            };

        mesh_run_nonloc(mesh, [&filler_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL)
                              {
                                  filler_nonloc({iL, node_info::X}, {jNL, node_info::X}, elL, elNL);
                                  filler_nonloc({iL, node_info::X}, {jNL, node_info::Y}, elL, elNL);
                                  filler_nonloc({iL, node_info::Y}, {jNL, node_info::X}, elL, elNL);
                                  filler_nonloc({iL, node_info::Y}, {jNL, node_info::Y}, elL, elNL);
                              });
    }

    return {std::move(triplets), std::move(triplets_bound)};
}

static void create_matrix(const mesh_2d<double> &mesh, const parameters<double> &params,
                          const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                                       boundary_type, std::function<double(double, double)>>> &bounds_cond,
                          Eigen::SparseMatrix<double> &K, Eigen::SparseMatrix<double> &K_bound,
                          const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    double time = omp_get_wtime();
    auto [triplets, triplets_bound] = triplets_fill(mesh, kinematic_nodes_set(mesh, bounds_cond), params, p1, influence_fun);
    std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;

    K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
    triplets_bound.clear();
    K.setFromTriplets(triplets.cbegin(), triplets.cend());
    std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
}

std::array<std::vector<double>, 3> strains_calc(const mesh_2d<double> &mesh, const Eigen::VectorXd &u)
{
    std::vector<double> eps11(mesh.nodes_count()),
                        eps22(mesh.nodes_count()),
                        eps12(mesh.nodes_count());

    std::array<double, 4> jacobi;
    std::vector<uint8_t> repeating(mesh.nodes_count(), 0);
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            ++repeating[mesh.node_number(el, i)];
            const std::array<double, 2> &node = e->node(i);
            memset(jacobi.data(), 0, jacobi.size() * sizeof(double));
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                jacobi[0] += mesh.coord(mesh.node_number(el, j), 0) * e->Nxi (j, node[0], node[1]);
                jacobi[1] += mesh.coord(mesh.node_number(el, j), 0) * e->Neta(j, node[0], node[1]);
                jacobi[2] += mesh.coord(mesh.node_number(el, j), 1) * e->Nxi (j, node[0], node[1]);
                jacobi[3] += mesh.coord(mesh.node_number(el, j), 1) * e->Neta(j, node[0], node[1]);
            }

            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                double jacobian = jacobi[0]*jacobi[3] - jacobi[1]*jacobi[2],
                       dx1 =  jacobi[3] * e->Nxi(j, node[0], node[1]) - jacobi[2] * e->Neta(j, node[0], node[1]),
                       dx2 = -jacobi[1] * e->Nxi(j, node[0], node[1]) + jacobi[0] * e->Neta(j, node[0], node[1]);
                eps11[mesh.node_number(el, i)] +=  dx1 * u[2*mesh.node_number(el, j)  ]  / jacobian;
                eps22[mesh.node_number(el, i)] +=  dx2 * u[2*mesh.node_number(el, j)+1]  / jacobian;
                eps12[mesh.node_number(el, i)] += (dx2 * u[2*mesh.node_number(el, j)  ] +
                                                   dx1 * u[2*mesh.node_number(el, j)+1]) / jacobian;
            }
        }
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        eps11[i] /=   repeating[i];
        eps22[i] /=   repeating[i];
        eps12[i] /= 2*repeating[i];
    }

    return {std::move(eps11), std::move(eps22), std::move(eps12)};
}

// Получение напряжений в узлах сетки, путём их переинтерполяции из квадратурных узлов
// Данный кусок взят из другой программы. Я не до конца понимаю как это работает, на мой взгляд это работать не должно.
// Пока что будем считать, что сетка однородная и состоит из билинейных элементов. В будущем, я надеюсь, это будет исправлено.
std::array<std::vector<double>, 3> stress_calc(const mesh_2d<double> &mesh, const Eigen::VectorXd &u, const parameters<double> &params)
{
    std::vector<double> sigma11(mesh.nodes_count()),
                        sigma22(mesh.nodes_count()),
                        sigma12(mesh.nodes_count());

    std::vector<uint8_t> repeating(mesh.nodes_count(), 0);
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ++repeating[mesh.node_number(el, i)];
    }

    Eigen::MatrixXd NQP(e->nodes_count(), e->nodes_count());
    for(size_t q = 0; q < e->nodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
            NQP(q, i) = e->qN(i, q);
    Eigen::MatrixXd NQPI = NQP.inverse();

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
    D(0, 0) = D(1, 1) = params.E / (1. - params.nu*params.nu);
    D(0, 1) = D(1, 0) = params.nu * params.E / (1. - params.nu*params.nu);
    D(2, 2) = 0.5 * params.E / (1. + params.nu);

    Eigen::VectorXd SigmaXX_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaYY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaXY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaXX_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaYY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),	
	                SigmaXY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                Ue              = Eigen::VectorXd::Zero(2 * e->nodes_count()),
	                Epsi            = Eigen::VectorXd::Zero(3),
                    Sigma           = Eigen::VectorXd::Zero(3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count()),
                    ElementNodesCoord(e->nodes_count(), 2),
                    Ndx(2, e->nodes_count());
    Eigen::Matrix2d Jmatr;
    Eigen::RowVectorXi ElementNodesNumbers(e->nodes_count());
    std::vector<Eigen::MatrixXd> NGradArr(e->nodes_count(), Eigen::MatrixXd::Zero(2, e->nodes_count()));

    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count());
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ElementNodesNumbers[i] = mesh.node_number(el, i);

        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            ElementNodesCoord(i, 0) = mesh.coord(ElementNodesNumbers[i], 0);
            ElementNodesCoord(i, 1) = mesh.coord(ElementNodesNumbers[i], 1);
            Ue(i * 2)     = u(ElementNodesNumbers(i) * 2);
			Ue(i * 2 + 1) = u(ElementNodesNumbers(i) * 2 + 1);
        }

        for(size_t q = 0; q < e->nodes_count(); ++q)
        {
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                NGradArr[q](0, j) = e->qNxi(j, q);
                NGradArr[q](1, j) = e->qNeta(j, q);
            }

            Jmatr = NGradArr[q] * ElementNodesCoord;
			Ndx = Jmatr.inverse() * NGradArr[q];
			for(size_t k = 0; k < e->nodes_count(); ++k)
			{
				B(0, k * 2)     = Ndx(0, k);
				B(1, k * 2 + 1) = Ndx(1, k);
				B(2, k * 2)     = Ndx(1, k);
                B(2, k * 2 + 1) = Ndx(0, k);
			}

            Epsi =  B * Ue;
			Epsi(2) *= 0.5;
			Sigma = D * Epsi;
			
			SigmaXX_QP(q) = Sigma(0);
			SigmaYY_QP(q) = Sigma(1);
			SigmaXY_QP(q) = Sigma(2);
        }

        SigmaXX_Element = NQPI * SigmaXX_QP;
		SigmaYY_Element = NQPI * SigmaYY_QP;
		SigmaXY_Element = NQPI * SigmaXY_QP;

        for(size_t i = 0; i < e->nodes_count(); ++i)
		{
			sigma11[ElementNodesNumbers(i)] += SigmaXX_Element(i);
			sigma22[ElementNodesNumbers(i)] += SigmaYY_Element(i);
			sigma12[ElementNodesNumbers(i)] += SigmaXY_Element(i);
		}
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        sigma11[i] /= repeating[i];
        sigma22[i] /= repeating[i];
        sigma12[i] /= repeating[i];
    }

    return {std::move(sigma11), std::move(sigma22), std::move(sigma12)};
}

void stationary(const std::string &path, const mesh_2d<double> &mesh, const parameters<double> &params,
                const std::vector<std::tuple<boundary_type, std::function<double(double, double)>,
                                             boundary_type, std::function<double(double, double)>>> &bounds_cond,
                const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    Eigen::VectorXd f = Eigen::VectorXd::Zero(2*mesh.nodes_count());;
    Eigen::SparseMatrix<double> K(2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    
    double time = omp_get_wtime();
    create_matrix(mesh, params, bounds_cond, K, K_bound, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    boundary_condition(mesh, kinematic_nodes_vectors(mesh, bounds_cond), bounds_cond, 1., K_bound, f);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    //Eigen::PardisoLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    solver.compute(K);
    const Eigen::VectorXd u = solver.solve(f);
    //
    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;
  
    const auto [eps11, eps22, eps12] = strains_calc(mesh, u);
    const auto [sigma11, sigma22, sigma12] = stress_calc(mesh, u, params);

    save_as_vtk(path, mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);

    // RAW OUTPUT
    std::ofstream fout_ux(std::string("results//u_x.csv")),
                  fout_uy(std::string("results//u_y.csv"));
    fout_ux.precision(20);
    fout_uy.precision(20);
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        fout_ux << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i) << std::endl;
        fout_uy << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << u(2*i+1) << std::endl;
    }

    std::ofstream fout_eps11(std::string("results//eps11.csv")),
                  fout_eps22(std::string("results//eps22.csv")),
                  fout_eps12(std::string("results//eps12.csv"));
    fout_eps11.precision(20);
    fout_eps22.precision(20);
    fout_eps12.precision(20);
    for(size_t i = 0; i < eps11.size(); ++i)
    {
        fout_eps11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps11[i] << std::endl;
        fout_eps12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps12[i] << std::endl;
        fout_eps22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << eps22[i] << std::endl;
    }

    std::ofstream fout_sigma11(std::string("results//sigma11.csv")),
                  fout_sigma22(std::string("results//sigma22.csv")),
                  fout_sigma12(std::string("results//sigma12.csv"));
    fout_sigma11.precision(20);
    fout_sigma22.precision(20);
    fout_sigma12.precision(20);
    for(size_t i = 0; i < sigma11.size(); ++i)
    {
        fout_sigma11 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma11[i] << std::endl;
        fout_sigma12 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma12[i] << std::endl;
        fout_sigma22 << mesh.coord(i, 0) << "," << mesh.coord(i, 1) << "," << sigma22[i] << std::endl;
    }
}

}