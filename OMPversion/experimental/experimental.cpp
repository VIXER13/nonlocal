// В данном модуле реализуются экспериментальные алгоритмы, которые возможно в будущем будут встроены в программу
#include <set>
#include <algorithm>
#include <cassert>
#include <cstring>

#include "omp.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "power.hpp"
#include "finite_element_routine.hpp"

template<class Type, class Index>
static Type integrate_taylor(const mesh_2d<Type, Index> &mesh,
                             const size_t elL, const size_t iL, const size_t jL,
                             const std::vector<uint32_t> &shifts, const matrix<Type> &coords, const matrix<Type> &jacobi_matrices,
                             const Type p1, const std::function<Type(Type, Type, Type, Type)> &influence_fun)
{
    Type integralL = 0.0;
    auto *eL = dynamic_cast<const finite_element::element_2d_integrate<Type, finite_element::qubic_serendip>*>(mesh.element_2d(mesh.element_type(elL)));
    for(size_t q = 0, shiftL = shifts[elL]; q < eL->qnodes_count(); ++q, ++shiftL)
        integralL += (( eL->qNxi(iL, q)*jacobi_matrices(shiftL, 3) - eL->qNeta(iL, q)*jacobi_matrices(shiftL, 2)) *
                      ( eL->qNxi(jL, q)*jacobi_matrices(shiftL, 3) - eL->qNeta(jL, q)*jacobi_matrices(shiftL, 2)) +
                      (-eL->qNxi(iL, q)*jacobi_matrices(shiftL, 1) + eL->qNeta(iL, q)*jacobi_matrices(shiftL, 0)) *
                      (-eL->qNxi(jL, q)*jacobi_matrices(shiftL, 1) + eL->qNeta(jL, q)*jacobi_matrices(shiftL, 0))) /
                     (jacobi_matrices(shiftL, 0)*jacobi_matrices(shiftL, 3) - jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 2)) * eL->weight(q);

    std::vector<Type> dTdx  (eL->qnodes_count()), dTdy    (eL->qnodes_count()),
                      d2Tdx2(eL->qnodes_count()), d2Tdxdy (eL->qnodes_count()), d2Tdy2  (eL->qnodes_count()),
                      d3Tdx3(eL->qnodes_count()), d3Tdx2dy(eL->qnodes_count()), d3Tdxdy2(eL->qnodes_count()), d3Tdy3(eL->qnodes_count());
    for(size_t qL = 0, shiftL = shifts[elL]; qL < eL->qnodes_count(); ++qL, ++shiftL)
    {
        Type jacobian = jacobi_matrices(shiftL, 0)*jacobi_matrices(shiftL, 3) - jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 2);

        dTdx    [qL] = ( eL->qNxi(jL, qL)*jacobi_matrices(shiftL, 3) - eL->qNeta(jL, qL)*jacobi_matrices(shiftL, 2)) / jacobian;
        dTdy    [qL] = (-eL->qNxi(jL, qL)*jacobi_matrices(shiftL, 1) + eL->qNeta(jL, qL)*jacobi_matrices(shiftL, 0)) / jacobian;

        d2Tdx2  [qL] = (     eL->qNxi2  (jL, qL) * jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 3) -
                        2. * eL->qNxieta(jL, qL) * jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 2) +
                             eL->qNeta2 (jL, qL) * jacobi_matrices(shiftL, 2)*jacobi_matrices(shiftL, 2)) / (jacobian * jacobian);

        d2Tdxdy [qL] = (    -eL->qNxi2  (jL, qL) * jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 1) +
                             eL->qNxieta(jL, qL) *(jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 0) + jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 2)) +
                            -eL->qNeta2 (jL, qL) * jacobi_matrices(shiftL, 2)*jacobi_matrices(shiftL, 0)) / (jacobian * jacobian);

        d2Tdy2  [qL] = (     eL->qNxi2  (jL, qL) * jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 1) -
                        2. * eL->qNxieta(jL, qL) * jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 0) +
                             eL->qNeta2 (jL, qL) * jacobi_matrices(shiftL, 0)*jacobi_matrices(shiftL, 0)) / (jacobian * jacobian);

        d3Tdx3  [qL] = (     eL->qNxi3   (jL, qL) * math_meta::power<3>(jacobi_matrices(shiftL, 3)) -
                        3. * eL->qNxi2eta(jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 3))*jacobi_matrices(shiftL, 2) +
                        3. * eL->qNxieta2(jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 2))*jacobi_matrices(shiftL, 3) -
                             eL->qNeta   (jL, qL) * math_meta::power<3>(jacobi_matrices(shiftL, 2))) / math_meta::power<3>(jacobian);

        d3Tdx2dy[qL] = (    -eL->qNxi3   (jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 3))*jacobi_matrices(shiftL, 1) +
                             eL->qNxi2eta(jL, qL) * jacobi_matrices(shiftL, 3) * (2.*jacobi_matrices(shiftL, 2)*jacobi_matrices(shiftL, 1) + jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 0)) -
                             eL->qNxieta2(jL, qL) * jacobi_matrices(shiftL, 2) * (2.*jacobi_matrices(shiftL, 3)*jacobi_matrices(shiftL, 0) + jacobi_matrices(shiftL, 2)*jacobi_matrices(shiftL, 1)) +
                             eL->qNeta   (jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 3))*jacobi_matrices(shiftL, 0)) / math_meta::power<3>(jacobian);

        d3Tdxdy2[qL] = (     eL->qNxi3   (jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 1))*jacobi_matrices(shiftL, 3) -
                             eL->qNxi2eta(jL, qL) * jacobi_matrices(shiftL, 1) * (2.*jacobi_matrices(shiftL, 0)*jacobi_matrices(shiftL, 3) + jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 2)) +
                             eL->qNxieta2(jL, qL) * jacobi_matrices(shiftL, 0) * (2.*jacobi_matrices(shiftL, 1)*jacobi_matrices(shiftL, 2) + jacobi_matrices(shiftL, 0)*jacobi_matrices(shiftL, 3)) -
                             eL->qNeta   (jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 0))*jacobi_matrices(shiftL, 2)) / math_meta::power<3>(jacobian);

        d3Tdy3  [qL] = (    -eL->qNxi3   (jL, qL) * math_meta::power<3>(jacobi_matrices(shiftL, 1)) +
                        3. * eL->qNxi2eta(jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 1))*jacobi_matrices(shiftL, 0) -
                        3. * eL->qNxieta2(jL, qL) * math_meta::power<2>(jacobi_matrices(shiftL, 0))*jacobi_matrices(shiftL, 1) +
                             eL->qNeta   (jL, qL) * math_meta::power<3>(jacobi_matrices(shiftL, 0))) / math_meta::power<3>(jacobian);
    }

    Type integralNL = 0.0,
         diff_x = 0., diff_y = 0., finit = 0.,
         int_const = 0., int_linear_x = 0., int_linear_y = 0., int_quad_xx = 0., int_quad_xy = 0., int_quad_yy = 0.;
    const finite_element::element_2d_integrate<Type, finite_element::qubic_serendip> *eNL = nullptr;
    for(auto elNL : mesh.neighbor(elL))
    {
        eNL = reinterpret_cast<const finite_element::element_2d_integrate<Type, finite_element::qubic_serendip>*>(mesh.element_2d(mesh.element_type(elNL)));
        for(size_t qL = 0, shiftL = shifts[elL]; qL < eL->qnodes_count(); ++qL, ++shiftL)
        {
            int_const = int_linear_x = int_linear_y = int_quad_xx = int_quad_xy = int_quad_yy = 0.;
            for(size_t qNL = 0, shiftNL = shifts[elNL]; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL)
            {
                diff_x = coords(shiftNL, 0) - coords(shiftL, 0),
                diff_y = coords(shiftNL, 1) - coords(shiftL, 1);
                finit = eNL->weight(qNL) * influence_fun(coords(shiftL, 0), coords(shiftNL, 0), coords(shiftL, 1), coords(shiftNL, 1)) * 
                        (jacobi_matrices(shiftNL, 0)*jacobi_matrices(shiftNL, 3) - jacobi_matrices(shiftNL, 1)*jacobi_matrices(shiftNL, 2));
                int_const    += finit;
                int_linear_x += finit * diff_x;
                int_linear_y += finit * diff_y;
                int_quad_xx  += finit * diff_x * diff_x;
                int_quad_xy  += finit * diff_x * diff_y;
                int_quad_yy  += finit * diff_y * diff_y;
            }
            integralNL += eL->weight(qL) * 
                          (( eL->qNxi(iL, qL)*jacobi_matrices(shiftL, 3) - eL->qNeta(iL, qL)*jacobi_matrices(shiftL, 2)) *
                           (dTdx[qL]*int_const /*+ d2Tdx2 [qL]*int_linear_x + d2Tdxdy[qL]*int_linear_y + 0.5*d3Tdx3  [qL]*int_quad_xx + d3Tdx2dy[qL]*int_quad_xy + 0.5*d3Tdxdy2[qL]*int_quad_yy*/) +
                           (-eL->qNxi(iL, qL)*jacobi_matrices(shiftL, 1) + eL->qNeta(iL, qL)*jacobi_matrices(shiftL, 0)) *
                           (dTdy[qL]*int_const /*+ d2Tdxdy[qL]*int_linear_x + d2Tdy2 [qL]*int_linear_y + 0.5*d3Tdx2dy[qL]*int_quad_xx + d3Tdxdy2[qL]*int_quad_xy + 0.5*d3Tdy3  [qL]*int_quad_yy)*/));
        }
    }

    return p1 * integralL + (1.-p1) * integralNL;
}

static void triplets_run_taylor(const mesh_2d<double> &mesh, std::vector<Eigen::Triplet<double>> &triplets, const size_t start,
                                const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    size_t el = 0;
    std::vector<uint32_t> shifts = quadrature_shifts_init(mesh);
    matrix<double> all_quad_coords = approx_all_quad_nodes_coords(mesh, shifts),
                   all_jacobi_matrices = approx_all_jacobi_matrices(mesh, shifts);
    for(size_t i = start; i < triplets.size(); ++i)
    {
        el = *reinterpret_cast<const size_t*>(&triplets[i].value());
        triplets[i] = Eigen::Triplet<double>(mesh.node_number(el, triplets[i].row()),
                                             mesh.node_number(el, triplets[i].col()),
                                             integrate_taylor<double>(mesh, el, triplets[i].row(), triplets[i].col(), 
                                                                      shifts, all_quad_coords, all_jacobi_matrices, p1, influence_fun));
    }
}

static void create_matrix_taylor(const mesh_2d<double> &mesh,
                                 const std::vector<std::tuple<boundary_type, std::function<double(double, double)>>> &bounds_cond,
                                 Eigen::SparseMatrix<double> &K, Eigen::SparseMatrix<double> &K_bound,
                                 const double p1, const std::function<double(double, double, double, double)> &influence_fun)
{
    std::set<uint32_t> temperature_nodes = temperature_nodes_set(mesh, bounds_cond);
    auto [triplets, classic_count, triplets_bound, classic_bound_count] = mesh_analysis(mesh, temperature_nodes, false);

    double time = omp_get_wtime();

    triplets_run_taylor(mesh, triplets, temperature_nodes.size(), p1, influence_fun);
    triplets_run_taylor(mesh, triplets_bound, 0, p1, influence_fun);
    
    std::cout << "Triplets calc: " << omp_get_wtime() - time << std::endl;

    temperature_nodes.clear();
    K_bound.setFromTriplets(triplets_bound.begin(), triplets_bound.end());
    triplets_bound.clear();
    K.setFromTriplets(triplets.begin(), triplets.end());
}