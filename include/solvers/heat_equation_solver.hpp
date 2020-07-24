#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <iostream>
#include <algorithm>
#include <omp.h>
#include "finite_element_routine.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/PardisoSupport"

namespace nonlocal::heat {

enum class boundary_t : uint8_t {
    TEMPERATURE = uint8_t(boundary_type::FIRST_KIND),
    FLOW        = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
struct boundary_condition {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    std::function<T(const std::array<T, 2>&)> func = [](const std::array<T, 2>&) noexcept { return 0; };
    boundary_t type = boundary_t::FLOW;
};

class _heat : protected _finite_element_routine {
protected:
    explicit _heat() noexcept = default;

    // Интегрирование базисной функции i по элементу e.
    // Для использования, предварительно должны быть проинициализированы jacobi_matrices текущего элемента.
    // Квадратурный сдвиг quad_shift нужен в случае если матрицы Якоби текуще элемента хранятся со сдвигом.
    template<class T, class Finite_Element_2D_Ptr>
    static T integrate_basic(const Finite_Element_2D_Ptr& e, const size_t i, 
                             const std::vector<std::array<T, 4>>& jacobi_matrices, size_t quad_shift) {
        T integral = 0;
        for(size_t q = 0; q < e->nodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) * e->qN(i, q) * jacobian(jacobi_matrices[quad_shift]);
        return integral;
    }

    // Интегрирование произведения базисных функций i и j элемента e.
    template<class T, class Finite_Element_2D_Ptr>
    static T integrate_basic_pair(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, 
                                  const std::vector<std::array<T, 4>>& jacobi_matrices, size_t quad_shift) {
        T integral = 0;
        for(size_t q = 0; q < e->nodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) * e->qN(i, q) * e->qN(j, q) * jacobian(jacobi_matrices[quad_shift]);
        return integral;
    }

    // Интегрирование произведения пар градиентов функций i и j элемента e.
    template<class T, class Finite_Element_2D_Ptr>
    static T integrate_gradient_pair(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, 
                                     const std::vector<std::array<T, 4>>& jacobi_matrices, size_t quad_shift) {
        T integral = 0;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += (( e->qNxi(i, q) * jacobi_matrices[quad_shift][3] - e->qNeta(i, q) * jacobi_matrices[quad_shift][2]) *
                         ( e->qNxi(j, q) * jacobi_matrices[quad_shift][3] - e->qNeta(j, q) * jacobi_matrices[quad_shift][2]) +
                         (-e->qNxi(i, q) * jacobi_matrices[quad_shift][1] + e->qNeta(i, q) * jacobi_matrices[quad_shift][0]) *
                         (-e->qNxi(j, q) * jacobi_matrices[quad_shift][1] + e->qNeta(j, q) * jacobi_matrices[quad_shift][0])) *
                        e->weight(q) / jacobian(jacobi_matrices[quad_shift]);
        return integral;
    }

    // Интегрирование произведения пары градиентов функции iL элемента eL и функции jNL элемента eNL.
    // Influence_Function - функтор с сигнатурой T(std::array<T, 2>&, std::array<T, 2>&)
    template<class T, class Finite_Element_2D_Ptr, class Influence_Function>
    static T integrate_gradient_pair_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,
                                            const Finite_Element_2D_Ptr& eNL, const size_t jNL,
                                            const std::vector<std::array<T, 2>>& quad_coords,
                                            const std::vector<std::array<T, 4>>& jacobi_matrices,
                                            size_t shiftL, size_t shiftNL, const Influence_Function& influence_function) {
        T integral = 0;
        const size_t sub_shift = shiftNL;
        for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL) {
            T int_with_weight_x = 0, int_with_weight_y = 0;
            for(size_t qNL = 0, shiftNL = sub_shift; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL) {
                const T influence_weight = eNL->weight(qNL) * influence_function(quad_coords[shiftL], quad_coords[shiftNL]);
                int_with_weight_x += influence_weight * ( eNL->qNxi (jNL, qNL) * jacobi_matrices[shiftNL][3] - 
                                                          eNL->qNeta(jNL, qNL) * jacobi_matrices[shiftNL][2]);
                int_with_weight_y += influence_weight * (-eNL->qNxi (jNL, qNL) * jacobi_matrices[shiftNL][1] + 
                                                          eNL->qNeta(jNL, qNL) * jacobi_matrices[shiftNL][0]);
            }
            integral += eL->weight(qL) *
                        (int_with_weight_x * ( eL->qNxi (iL, qL) * jacobi_matrices[shiftL][3] -
                                               eL->qNeta(iL, qL) * jacobi_matrices[shiftL][2]) +
                         int_with_weight_y * (-eL->qNxi (iL, qL) * jacobi_matrices[shiftL][1] + 
                                               eL->qNeta(iL, qL) * jacobi_matrices[shiftL][0]));
        }
        return integral;
    }

    template<class Type, class Index, class Vector, class Right_Part>
    static void integrate_right_part(Vector& f, const mesh::mesh_2d<Type, Index>& mesh, 
                                     const std::vector<Index>& shifts_quad,
                                     const std::vector<std::array<Type, 2>>& all_quad_coords,
                                     const std::vector<std::array<Type, 4>>& all_jacobi_matrices,
                                     const Right_Part& right_part) {
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t i = 0; i < e->nodes_count(); ++i)
                f[mesh.node_number(el, i)] += integrate_function(e, i, all_quad_coords, all_jacobi_matrices, shifts_quad[el], right_part);
        }
    }

    // Функция возвращает массив флагов, где true свидетельствует о том, что узел под данным номером внутренний,
    // т.е. на нём не задано граничное условие первого рода
    template<class Type, class Index>
    static std::vector<bool> inner_nodes_vector(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond) {
        std::vector<bool> inner_nodes(mesh.nodes_count(), true);
        for(size_t b = 0; b < mesh.boundary_groups_count(); ++b)
            if(bounds_cond[b].type == boundary_t::TEMPERATURE)
                for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                    const auto& e = mesh.element_1d(mesh.element_1d_type(b, el));
                    for(size_t i = 0; i < e->nodes_count(); ++i)
                        inner_nodes[mesh.node_number(b, el, i)] = false;
                }
        return std::move(inner_nodes);
    }

    // Создаёт массив векторов размерности равной количеству граничных условий и записывает в них узлы,
    // на которых заданы граничные условия первого рода, исключая повторяющиеся.
    template<class Type, class Index>
    static std::vector<std::vector<Index>> temperature_nodes_vectors(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond) {
        std::vector<std::vector<Index>> temperature_nodes(mesh.boundary_groups_count());
        for(size_t b = 0; b < mesh.boundary_groups_count(); ++b)
            if(bounds_cond[b].type == boundary_t::TEMPERATURE)
                for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                    const auto& e = mesh.element_1d(mesh.element_1d_type(b, el));
                    for(size_t i = 0; i < e->nodes_count(); ++i) {
                        bool push = true;
                        for(const auto& bound : temperature_nodes)
                            push = push && std::find(bound.cbegin(), bound.cend(), mesh.node_number(b, el, i)) == bound.cend();
                        if(push)
                            temperature_nodes[b].push_back(mesh.node_number(b, el, i));
                    }
                }
        return std::move(temperature_nodes);
    }

    // Данная функция анализирует сетку и вычисляет сдвиги для дальнейшего интегрирования.
    // Здесь же происходит расщепление итоговой матрицы теплопроводности на две части:
    // первая, которая будет представлять из себя СЛАУ, а вторая, которая перейдёт в правую часть.
    template<class Type, class Index>
    static std::array<std::vector<Index>, 4>
        mesh_analysis(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<bool>& inner_nodes, const bool nonlocal) {
        std::vector<Index> shifts_loc      (mesh.elements_count()+1, 0),
                           shifts_bound_loc(mesh.elements_count()+1, 0),
                           shifts_nonloc, shifts_bound_nonloc;

        mesh_run_loc(mesh,
            [&mesh, &inner_nodes, &shifts_loc, &shifts_bound_loc](size_t i, size_t j, size_t el) { 
                if(mesh.node_number(el, i) >= mesh.node_number(el, j)) {
                    if(inner_nodes[mesh.node_number(el, i)] && inner_nodes[mesh.node_number(el, j)])
                        ++shifts_loc[el+1];
                    else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                        ++shifts_bound_loc[el+1];
                }
            });

        shifts_loc[0] = std::count(inner_nodes.cbegin(), inner_nodes.cend(), false);
        for(size_t i = 1; i < mesh.elements_count()+1; ++i) {
            shifts_loc[i] += shifts_loc[i-1];
            shifts_bound_loc[i] += shifts_bound_loc[i-1];
        }
        
        if(nonlocal) {
            shifts_nonloc.resize(mesh.elements_count()+1, 0);
            shifts_bound_nonloc.resize(mesh.elements_count()+1, 0);

            mesh_run_nonloc(mesh, 
                [&mesh, &inner_nodes, &shifts_nonloc, &shifts_bound_nonloc](size_t iL, size_t jNL, size_t elL, size_t elNL) { 
                    if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL)) {
                        if(inner_nodes[mesh.node_number(elL, iL)] && inner_nodes[mesh.node_number(elNL, jNL)])
                            ++shifts_nonloc[elL+1];
                        else if(mesh.node_number(elL, iL) != mesh.node_number(elNL, jNL))
                            ++shifts_bound_nonloc[elL+1];
                    }
                });

            shifts_nonloc[0] = shifts_loc.back();
            shifts_bound_nonloc[0] = shifts_bound_loc.back();
            for(size_t i = 1; i < mesh.elements_count()+1; ++i) {
                shifts_nonloc[i] += shifts_nonloc[i-1];
                shifts_bound_nonloc[i] += shifts_bound_nonloc[i-1];
            }
        }

        return {std::move(shifts_loc), std::move(shifts_bound_loc), std::move(shifts_nonloc), std::move(shifts_bound_nonloc)};
    }

    // Функция заполнения триплетов, перед их сборкой в итоговую матрицу.
    // Integrate_Rule - функтор с сигнатурой Type(const Finite_Element_2D_Ptr&, const size_t, const size_t, 
    //                                            const std::vector<std::array<T, 4>>&, size_t)
    // Influence_Function - функтор с сигнатурой T(std::array<T, 2>&, std::array<T, 2>&)
    template<class Type, class Index, class Integrate_Rule, class Influence_Function>
    static std::array<std::vector<Eigen::Triplet<Type, Index>>, 2>
        triplets_fill(const mesh::mesh_2d<Type, Index>& mesh, 
                      const std::vector<Index>& shifts_quad,
                      const std::vector<std::array<Type, 2>>& all_quad_coords,
                      const std::vector<std::array<Type, 4>>& all_jacobi_matrices,
                      const std::vector<boundary_condition<Type>>& bounds_cond, const bool neumann_task,
                      const Integrate_Rule& integrate_rule,
                      const Type p1, const Influence_Function& influence_fun) {
        static constexpr Type MAX_LOCAL_WEIGHT = 0.999;
        const bool nonlocal = p1 < MAX_LOCAL_WEIGHT;
        const std::vector<bool> inner_nodes = inner_nodes_vector(mesh, bounds_cond);
        auto [shifts_loc, shifts_bound_loc, shifts_nonloc, shifts_bound_nonloc] = mesh_analysis(mesh, inner_nodes, nonlocal);

        size_t neumann_triplets = 0;
        if(neumann_task) {
            for(size_t el = 0; el < mesh.elements_count(); ++el)
                neumann_triplets += mesh.element_2d(mesh.element_2d_type(el))->nodes_count();
        }

        std::vector<Eigen::Triplet<Type, Index>> triplets      ((nonlocal ? shifts_nonloc.back()       : shifts_loc.back()) + neumann_triplets),
                                                 triplets_bound( nonlocal ? shifts_bound_nonloc.back() : shifts_bound_loc.back());
        if(!neumann_task)
            for(size_t i = 0, j = 0; i < inner_nodes.size(); ++i)
                if(!inner_nodes[i])
                    triplets[j++] = Eigen::Triplet<Type, Index>(i, i, 1.);

        mesh_run_loc(mesh,
            [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_loc, &shifts_bound_loc,
            &shifts_quad, &all_jacobi_matrices, &integrate_rule, p1] (size_t i, size_t j, size_t el) {
                if(mesh.node_number(el, i) >= mesh.node_number(el, j)) {
                    const Type integral = p1 * integrate_rule(mesh.element_2d(mesh.element_2d_type(el)), i, j, all_jacobi_matrices, shifts_quad[el]);
                    if(inner_nodes[mesh.node_number(el, i)] && inner_nodes[mesh.node_number(el, j)])
                        triplets[shifts_loc[el]++] = Eigen::Triplet<Type, Index>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
                    else if(mesh.node_number(el, i) != mesh.node_number(el, j))
                        triplets_bound[shifts_bound_loc[el]++] = inner_nodes[mesh.node_number(el, j)] ?
                                                                 Eigen::Triplet<Type, Index>(mesh.node_number(el, j), mesh.node_number(el, i), integral) :
                                                                 Eigen::Triplet<Type, Index>(mesh.node_number(el, i), mesh.node_number(el, j), integral);
                }
            });

        if(nonlocal) {
            mesh_run_nonloc(mesh, 
                [&mesh, &inner_nodes, &triplets, &triplets_bound, &shifts_nonloc, &shifts_bound_nonloc,
                &shifts_quad, &all_jacobi_matrices, &all_quad_coords, &influence_fun, p2 = 1. - p1]
                (size_t iL, size_t jNL, size_t elL, size_t elNL) {
                    if(mesh.node_number(elL, iL) >= mesh.node_number(elNL, jNL)) {
                        const Type integral = p2 * integrate_gradient_pair_nonloc(mesh.element_2d(mesh.element_2d_type(elL )), iL,
                                                                                  mesh.element_2d(mesh.element_2d_type(elNL)), jNL,
                                                                                  all_quad_coords, all_jacobi_matrices,
                                                                                  shifts_quad[elL], shifts_quad[elNL], influence_fun);
                        if(inner_nodes[mesh.node_number(elL, iL)] && inner_nodes[mesh.node_number(elNL, jNL)])
                            triplets[shifts_nonloc[elL]++] = Eigen::Triplet<Type, Index>(mesh.node_number(elL, iL), mesh.node_number(elNL, jNL), integral);
                        else if(mesh.node_number(elL, iL) != mesh.node_number(elNL, jNL))
                            triplets_bound[shifts_bound_nonloc[elL]++] = inner_nodes[mesh.node_number(elNL, jNL)] ?
                                                                         Eigen::Triplet<Type, Index>(mesh.node_number(elNL, jNL), mesh.node_number(elL,   iL), integral) :
                                                                         Eigen::Triplet<Type, Index>(mesh.node_number(elL,  iL ), mesh.node_number(elNL, jNL), integral);
                    }
                });
        }

        if(neumann_task)
        {
            size_t last_index = nonlocal ? shifts_nonloc.back() : shifts_loc.back();
            for(size_t el = 0; el < mesh.elements_count(); ++el) {
                const auto& e = mesh.element_2d(mesh.element_2d_type(el));
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    triplets[last_index++] = Eigen::Triplet<Type, Index>(mesh.nodes_count(), mesh.node_number(el, i), 
                                                                         integrate_basic(e, i, all_jacobi_matrices, shifts_quad[el]));
            }
        }

        return {std::move(triplets), std::move(triplets_bound)};
    }

    // Вычисление марицы теплопроводности (теплоёмкости в случае когда integrate_rule = integrate_basic_pair).
    // На выходе получаем расщеплённую матрицу, где K будет участвовать в решение СЛАУ, а K_bound уйдёт в правую часть.
    // Integrate_Rule - функтор с сигнатурой Type(const Finite_Element_2D_Ptr&, const size_t, const size_t, 
    //                                            const std::vector<std::array<T, 4>>&, size_t)
    // Influence_Function - функтор с сигнатурой Type(std::array<Type, 2>&, std::array<Type, 2>&)
    // P.S. На момент написания, в Eigen были проблемы с move-семантикой, поэтому вопреки выше описанным функциям,
    // матрицы K и K_bound передаются по ссылке в функцию.
    template<class Type, class Index, class Integrate_Rule, class Influence_Function>
    static void create_matrix(Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K, Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K_bound,
                              const mesh::mesh_2d<Type, Index>& mesh, 
                              const std::vector<Index>& shifts_quad,
                              const std::vector<std::array<Type, 2>>& all_quad_coords,
                              const std::vector<std::array<Type, 4>>& all_jacobi_matrices,
                              const std::vector<boundary_condition<Type>>& bounds_cond, const bool neumann_task,
                              const Integrate_Rule& integrate_rule, const Type p1, const Influence_Function& influence_fun) {
        double time = omp_get_wtime();
        std::cout << "Triplets calc: ";
        auto [triplets, triplets_bound] = triplets_fill(mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, 
                                                        bounds_cond, neumann_task, integrate_rule, p1, influence_fun);
        std::cout << omp_get_wtime() - time << std::endl
                  << "Triplets count: " << triplets.size() + triplets_bound.size() << std::endl;

        K_bound.setFromTriplets(triplets_bound.cbegin(), triplets_bound.cend());
        triplets_bound.reserve(0);
        K.setFromTriplets(triplets.cbegin(), triplets.cend());
        std::cout << "Nonzero elemets count: " << K.nonZeros() + K_bound.nonZeros() << std::endl;
    }

    template<class Type, class Index, class Vector>
    static void integrate_boundary_flow(Vector& f, const mesh::mesh_2d<Type, Index>& mesh,
                                        const std::vector<boundary_condition<Type>>& bounds_cond) {
        std::vector<std::array<Type, 2>> quad_nodes, jacobi_matrices;
        for(size_t b = 0; b < bounds_cond.size(); ++b)
            if(bounds_cond[b].type == boundary_t::FLOW)
                for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                    approx_quad_nodes_on_bound(quad_nodes, mesh, b, el);
                    approx_jacobi_matrices_on_bound(jacobi_matrices, mesh, b, el);
                    const auto& be = mesh.element_1d(mesh.element_1d_type(b, el));
                    for(size_t i = 0; i < be->nodes_count(); ++i)
                        f[mesh.node_number(b, el, i)] += 
                            integrate_boundary_gradient(be, i, quad_nodes, jacobi_matrices, bounds_cond[b].func);
                }
    }

    // Учёт граничных условий.
    // Параметр tau - шаг интегрирования, в стационарной задаче равен 1.
    template<class Type, class Index, class Vector>
    static void temperature_on_boundary(Vector& f, const mesh::mesh_2d<Type, Index>& mesh,
                                        const std::vector<boundary_condition<Type>>& bounds_cond,
                                        const Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>& K_bound) {
        // Граничные условия первого рода
        const std::vector<std::vector<Index>> temperature_nodes = temperature_nodes_vectors(mesh, bounds_cond);
        for(size_t b = 0; b < temperature_nodes.size(); ++b)
            for(const Index node : temperature_nodes[b]) {
                const Type temp = bounds_cond[b].func(mesh.node(node));
                for(typename Eigen::SparseMatrix<Type>::InnerIterator it(K_bound, node); it; ++it)
                    f[it.row()] -= temp * it.value();
            }

        // Повторный проход для корректировки
        for(size_t b = 0; b < temperature_nodes.size(); ++b)
            for(const Index node : temperature_nodes[b])
                f[node] = bounds_cond[b].func(mesh.node(node));
    }

public:
    template<class Type, class Index, class Vector>
    friend Type integrate_solution(const mesh::mesh_2d<Type, Index>& mesh, const Vector& T);

    // Функция, решающая стационарное уравнение теплопроводности в нелокальной постановке.
    // Right_Part - функтор с сигнатурой Type(std::array<Type, 2>&),
    // Influence_Function - функтор с сигнатурой Type(std::array<Type, 2>&, std::array<Type, 2>&)
    // volume - значение интеграла по области, в случае если поставлена задача Неймана, по умолчанию 0.
    template<class Type, class Index, class Right_Part, class Influence_Function>
    friend Eigen::Matrix<Type, Eigen::Dynamic, 1>
        stationary(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond, 
                   const Right_Part& right_part, const Type p1, const Influence_Function& influence_fun, const Type volume = 0);

    template<class Type, class Index, class Init_Distribution, class Right_Part, class Influence_Function>
    friend void nonstationary(const std::string& path, 
                              const mesh::mesh_2d<Type, Index>& mesh, const Type tau, const uintmax_t time_steps,
                              const std::vector<boundary_condition<Type>>& bounds_cond,
                              const Init_Distribution& init_dist,
                              const Right_Part& right_part,
                              const Type p1, const Influence_Function& influence_fun,
                              const uintmax_t print_frequency);
};

template<class Type, class Index, class Vector>
Type integrate_solution(const mesh::mesh_2d<Type, Index>& mesh, const Vector& T) {
    if(mesh.nodes_count() != size_t(T.size()))
        throw std::logic_error{"mesh.nodes_count() != T.size()"};
    Type integral = 0;
    const std::vector<Index> shifts_quad = _heat::quadrature_shifts_init(mesh);
    const std::vector<std::array<Type, 4>> all_jacobi_matrices = _heat::approx_all_jacobi_matrices(mesh, shifts_quad);
    for(size_t el = 0; el < mesh.elements_count(); ++el) {
        const auto& e = mesh.element_2d(mesh.element_2d_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t q = 0, shift = shifts_quad[el]; q < e->qnodes_count(); ++q, ++shift)
                integral += e->weight(q) * e->qN(i, q) * T[mesh.node_number(el, i)] * _heat::jacobian(all_jacobi_matrices[shift]);
    }
    return integral;
}

template<class Type, class Index, class Right_Part, class Influence_Function>
Eigen::Matrix<Type, Eigen::Dynamic, 1>
    stationary(const mesh::mesh_2d<Type, Index>& mesh, const std::vector<boundary_condition<Type>>& bounds_cond,
               const Right_Part& right_part, const Type p1, const Influence_Function& influence_fun, const Type volume) {
    const bool neumann_task = std::all_of(bounds_cond.cbegin(), bounds_cond.cend(),
                              [](const boundary_condition<Type>& bound) { return bound.type == boundary_t::FLOW; });
    const size_t matrix_size = neumann_task ? mesh.nodes_count()+1 : mesh.nodes_count();
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (matrix_size, matrix_size),
                                                      K_bound(matrix_size, matrix_size);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(matrix_size);

    _heat::integrate_boundary_flow<Type, Index>(f, mesh, bounds_cond);

    if(neumann_task) {
        Type sum = 0;
        for(size_t i = 0; i < matrix_size; ++i)
            sum += f[i];
        if(std::abs(sum) > 1e-5)
            throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
        f[mesh.nodes_count()] = volume;
    }

    std::cout << "Quadratures data init: ";
    double time = omp_get_wtime();
    const std::vector<Index> shifts_quad = _heat::quadrature_shifts_init(mesh);
    const std::vector<std::array<Type, 2>> all_quad_coords = _heat::approx_all_quad_nodes(mesh, shifts_quad);
    const std::vector<std::array<Type, 4>> all_jacobi_matrices = _heat::approx_all_jacobi_matrices(mesh, shifts_quad);
    std::cout << omp_get_wtime() - time << std::endl;

    std::cout << "Right part Integrate: ";
    time = omp_get_wtime();
    _heat::integrate_right_part(f, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, right_part);
    std::cout << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    _heat::create_matrix<Type, Index>(
        K, K_bound, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, bounds_cond, neumann_task,
        _heat::integrate_gradient_pair<Type, typename mesh::mesh_2d<Type, Index>::fe_2d_ptr>, 
        p1, influence_fun
    );
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    std::cout << "Boundary filling: ";
    time = omp_get_wtime();
    _heat::temperature_on_boundary<Type, Index>(f, mesh, bounds_cond, K_bound);
    std::cout << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> T = solver.solve(f);
    std::cout << "System solving: " << omp_get_wtime() - time << std::endl;
    T.conservativeResize(mesh.nodes_count());

    return std::move(T);
}

template<class Type, class Index, class Init_Distribution, class Right_Part, class Influence_Function>
void nonstationary(const std::string& path, 
                   const mesh::mesh_2d<Type, Index>& mesh, const Type tau, const uintmax_t time_steps,
                   const std::vector<boundary_condition<Type>>& bounds_cond,
                   const Init_Distribution& init_dist,
                   const Right_Part& right_part,
                   const Type p1, const Influence_Function& influence_fun,
                   const uintmax_t print_frequency) {
    const std::vector<Index> shifts_quad = _heat::quadrature_shifts_init(mesh);
    const std::vector<std::array<Type, 2>> all_quad_coords = _heat::approx_all_quad_nodes(mesh, shifts_quad);
    const std::vector<std::array<Type, 4>> all_jacobi_matrices = _heat::approx_all_jacobi_matrices(mesh, shifts_quad);
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (mesh.nodes_count(), mesh.nodes_count()),
                                                      K_bound(mesh.nodes_count(), mesh.nodes_count()),
                                                      C      (mesh.nodes_count(), mesh.nodes_count()),
                                                      C_bound(mesh.nodes_count(), mesh.nodes_count());

    static constexpr bool NOT_NEUMANN_TASK = false;
    _heat::create_matrix<Type, Index>(
        K, K_bound, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, bounds_cond, NOT_NEUMANN_TASK,
        _heat::integrate_gradient_pair<Type, typename mesh::mesh_2d<Type, Index>::fe_2d_ptr>, 
        p1, influence_fun
    );

    static constexpr Type LOCAL = 1;
    _heat::create_matrix<Type, Index>(
        C, C_bound, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, bounds_cond, NOT_NEUMANN_TASK,
        _heat::integrate_basic_pair<Type, typename mesh::mesh_2d<Type, Index>::fe_2d_ptr>, 
        LOCAL, influence_fun
    );

    C_bound.setZero();
    K_bound *= tau;
    K *= tau;
    K += C;
    for(size_t b = 0; b < mesh.boundary_groups_count(); ++b)
        if(bounds_cond[b].type == boundary_t::TEMPERATURE)
            for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                const auto& e = mesh.element_1d(mesh.element_1d_type(b, el));
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    K.coeffRef(mesh.node_number(b, el, i), mesh.node_number(b, el, i)) = 1;
            }

    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(mesh.nodes_count()),
                                           T_prev(mesh.nodes_count()), T(mesh.nodes_count());
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        T_prev[i] = init_dist(mesh.node(i));

    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);
    if(print_frequency != uintmax_t(-1)) {
        mesh.save_as_vtk(path + "0.vtk", T_prev);
        std::cout << "step = " << 0 << " Volume = " << integrate_solution(mesh, T_prev) << std::endl;
    }
    for(size_t i = 1; i < time_steps; ++i) {
        f.setZero();
        _heat::integrate_boundary_flow<Type, Index>(f, mesh, bounds_cond);
        _heat::integrate_right_part(f, mesh, shifts_quad, all_quad_coords, all_jacobi_matrices, right_part);
        f *= tau;
        f += C.template selfadjointView<Eigen::Lower>() * T_prev;
        _heat::temperature_on_boundary<Type, Index>(f, mesh, bounds_cond, K_bound);
        T = solver.solve(f);
        T_prev.swap(T);
        if(i % print_frequency == 0) {
            mesh.save_as_vtk(path + std::to_string(i) + ".vtk", T_prev);
            std::cout << "step = " << i << " Volume = " << integrate_solution(mesh, T_prev) << std::endl;
        }
    }
}

}

#endif