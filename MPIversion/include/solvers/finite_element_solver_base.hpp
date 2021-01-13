#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <petsc.h>
#include <petscsystypes.h>
#include "../../Eigen/Eigen/Sparse"
#include "mesh.hpp"
#include "utils.hpp"

namespace nonlocal {

enum class boundary_type : uint8_t {
    FIRST_KIND,
    SECOND_KIND
};

template<class T, class B, size_t N>
struct boundary_condition {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

    struct boundary_pair {
        B type = B(boundary_type::SECOND_KIND);
        std::function<T(const std::array<T, 2>&)> func = [](const std::array<T, 2>&) constexpr noexcept { return 0; };
    };

    std::array<boundary_pair, N> data;

    static constexpr size_t degrees_of_freedom() { return N; }

    B type(const size_t b) const { return data[b].type; }
    const std::function<T(const std::array<T, 2>&)>& func(const size_t b) const { return data[b].func; }

    bool contains_condition_second_kind() const {
        return std::any_of(data.cbegin(), data.cend(), [](const boundary_pair& pair) { return pair.type == B(boundary_type::SECOND_KIND); });
    }
};

template<class T, class I>
class finite_element_solver_base {
    mesh::mesh_2d<T, I>           _mesh;
    std::vector<I>                _quad_shifts;        // Квадратурные сдвиги
    std::vector<std::array<T, 2>> _quad_coords;        // Координаты квадратурных узлов сетки
    std::vector<std::array<T, 4>> _jacobi_matrices;    // Матрицы Якоби вычисленные в квадратурных узлах
    std::vector<std::vector<I>>   _nodes_elements_map; // Номера элементов, в которых присутствует узел
    std::vector<std::unordered_map<I, uint8_t>> _global_to_local_numbering; // Переход от глобальной нумерации к локальной каждого элемента
                                                                            // Считаем, что в элементе не более 255 узлов.
    std::vector<std::vector<I>>   _elements_neighbors; // Массив с номерами ближайших соседей

    // Квадратурные сдвиги по элементам.
    static std::vector<I> quadrature_shifts_init(const mesh::mesh_2d<T, I>& mesh) {
        std::vector<I> quad_shifts(mesh.elements_count()+1);
        quad_shifts[0] = 0;
        for(size_t el = 0; el < mesh.elements_count(); ++el)
            quad_shifts[el+1] = quad_shifts[el] + mesh.element_2d(mesh.element_2d_type(el))->qnodes_count();
        return std::move(quad_shifts);
    }

    // Аппроксимация глобальных координат всех квадратурных узлов сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    static std::vector<std::array<T, 2>> approx_all_quad_nodes(const mesh::mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts) {
        if(mesh.elements_count()+1 != quad_shifts.size())
            throw std::logic_error{"Quadrature shifts are incorrect."};
        std::vector<std::array<T, 2>> quad_coords(quad_shifts.back(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, quad_coords)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 2; ++comp)
                        quad_coords[quad_shifts[el]+q][comp] += mesh.node(mesh.node_number(el, i))[comp] * e->qN(i, q);
        }
        return std::move(quad_coords);
    }

    // Аппроксимация матриц Якоби во всех квадратурных узлах сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    static std::vector<std::array<T, 4>> approx_all_jacobi_matrices(const mesh::mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts) {
        if(mesh.elements_count()+1 != quad_shifts.size())
            throw std::logic_error{"Quadrature shifts are incorrect."};
        std::vector<std::array<T, 4>> jacobi_matrices(quad_shifts.back(), std::array<T, 4>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, jacobi_matrices)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                for(size_t i = 0; i < e->nodes_count(); ++i) {
                    jacobi_matrices[quad_shifts[el]+q][0] += mesh.node(mesh.node_number(el, i))[0] * e->qNxi (i, q);
                    jacobi_matrices[quad_shifts[el]+q][1] += mesh.node(mesh.node_number(el, i))[0] * e->qNeta(i, q);
                    jacobi_matrices[quad_shifts[el]+q][2] += mesh.node(mesh.node_number(el, i))[1] * e->qNxi (i, q);
                    jacobi_matrices[quad_shifts[el]+q][3] += mesh.node(mesh.node_number(el, i))[1] * e->qNeta(i, q);
                }
        }
        return std::move(jacobi_matrices);
    }

    static std::vector<std::vector<I>> node_elements_map_init(const mesh::mesh_2d<T, I>& mesh) {
        std::vector<std::vector<I>> nodes_elements_map(mesh.nodes_count());
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(el);
            for(size_t node = 0; node < e->nodes_count(); ++node)
                nodes_elements_map[mesh.node_number(el, node)].push_back(el);
        }
        return std::move(nodes_elements_map);
    }

    static std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering_init(const mesh::mesh_2d<T, I>& mesh) {
        std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh.elements_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(el);
            for(size_t node = 0; node < e->nodes_count(); ++node)
                global_to_local_numbering[el][mesh.node_number(el, node)] = node;
        }
        return std::move(global_to_local_numbering);
    }

    static std::vector<std::array<T, 2>> approx_centres_of_elements(const mesh::mesh_2d<T, I>& mesh) {
        std::vector<std::array<T, 2>> centres(mesh.elements_count(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, centres)
        for(size_t el = 0; el < centres.size(); ++el) {
            const auto& e = mesh.element_2d(el);
            const T x0 = mesh::is_trinagle(mesh.element_2d_type(el)) ? 1./3. : 0.;
            for(size_t node = 0; node < e->nodes_count(); ++node) {
                using namespace utils;
                centres[el] += mesh.node(mesh.node_number(el, node)) * e->N(node, {x0, x0});
            }
        }
        return std::move(centres);
    }

    static std::vector<std::vector<I>>
    find_elements_neighbors(const mesh::mesh_2d<T, I>& mesh, const std::vector<std::array<T, 2>>& centres,
                            const std::unordered_set<I>& elements, const T r) {
        std::vector<std::vector<I>> elements_neighbors(mesh.elements_count());
        for(const I elL : elements) {
            elements_neighbors[elL].reserve(mesh.elements_count());
            for(size_t elNL = 0; elNL < mesh.elements_count(); ++elNL)
                if(utils::distance(centres[elL], centres[elNL]) < r)
                    elements_neighbors[elL].push_back(elNL);
            elements_neighbors[elL].shrink_to_fit();
        }
        return std::move(elements_neighbors);
    }

    void init_first_and_last_nodes() {
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        _first_node = _mesh.nodes_count() / _size *  _rank;
        _last_node  = _mesh.nodes_count() / _size * (_rank+1) + (_rank == _size-1) * _mesh.nodes_count() % _size;
    }

protected:
    using Finite_Element_1D_Ptr = typename mesh::mesh_2d<T, I>::Finite_Element_1D_Ptr;
    using Finite_Element_2D_Ptr = typename mesh::mesh_2d<T, I>::Finite_Element_2D_Ptr;

    enum component : bool {X, Y};
    static constexpr T MAX_LOCAL_WEIGHT = 0.999;

    PetscMPIInt _size = 0, _rank = 0;
    size_t _first_node = 0, _last_node = 0;

    explicit finite_element_solver_base(const mesh::mesh_2d<T, I>& mesh) :
        _mesh{mesh},
        _quad_shifts{quadrature_shifts_init(_mesh)},
        _quad_coords{approx_all_quad_nodes(_mesh, _quad_shifts)},
        _jacobi_matrices{approx_all_jacobi_matrices(_mesh, _quad_shifts)},
        _nodes_elements_map{node_elements_map_init(_mesh)},
        _global_to_local_numbering{global_to_local_numbering_init(_mesh)} {
        init_first_and_last_nodes();

//        PetscMPIInt rank = -1, size = -1;
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        for(PetscMPIInt i = 0; i < size; ++i) {
//            if (i == rank)
//                std::cout << "_first_node = " << _first_node << std::endl
//                          << "_last_node  = " << _last_node  << std::endl;
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
    }

    explicit finite_element_solver_base(mesh::mesh_2d<T, I>&& mesh) :
        _mesh{std::move(mesh)},
        _quad_shifts{quadrature_shifts_init(_mesh)},
        _quad_coords{approx_all_quad_nodes(_mesh, _quad_shifts)},
        _jacobi_matrices{approx_all_jacobi_matrices(_mesh, _quad_shifts)},
        _nodes_elements_map{node_elements_map_init(_mesh)},
        _global_to_local_numbering{global_to_local_numbering_init(_mesh)} {
        init_first_and_last_nodes();
    }

    virtual ~finite_element_solver_base() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const { return _mesh; }
    I quad_shift(const size_t element) const { return _quad_shifts[element]; }
    const std::array<T, 2>& quad_coord(const size_t global_quad_node) const { return _quad_coords[global_quad_node]; }
    const std::array<T, 4>& jacobi_matrix(const size_t global_quad_node) const { return _jacobi_matrices[global_quad_node]; }
    const std::vector<I>& neighbors(const size_t element) const { return _elements_neighbors[element]; }
    const std::unordered_map<I, uint8_t>& global_to_local_numbering(const size_t element) const { return _global_to_local_numbering[element]; }
    const std::vector<I>& nodes_elements_map(const size_t node) const { return _nodes_elements_map[node]; }

    // Функция обхода сетки в локальных постановках.
    // Нужна для предварительного подсчёта количества элементов  и интегрирования системы.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t)
    template<class Callback>
    void mesh_run_loc(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback)
        for(size_t node = _first_node; node < _last_node; ++node) {
            for(const I el : _nodes_elements_map[node]) {
                const auto& e = _mesh.element_2d(el);
                const size_t i = _global_to_local_numbering[el].find(node)->second; // Проекционные функции
                for(size_t j = 0; j < e->nodes_count(); ++j) // Аппроксимационные функции
                    callback(el, i, j);
            }
        }
    }

    // Функция обхода сетки в нелокальных постановках.
    // Нужна для предварительного подсчёта количества элементов и интегрирования системы.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t, size_t)
    template<class Callback>
    void mesh_run_nonloc(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback)
        for(size_t node = _first_node; node < _last_node; ++node) {
            for(const I elL : _nodes_elements_map[node]) {
                const size_t iL =  _global_to_local_numbering[elL].find(node)->second; // Проекционные функции
                for(const I elNL : _elements_neighbors[elL]) {
                    const auto& eNL = _mesh.element_2d(elNL);
                    for(size_t jNL = 0; jNL < eNL->nodes_count(); ++jNL) // Аппроксимационные функции
                        callback(elL, iL, elNL, jNL);
                }
            }
        }
    }

    template<class Callback>
    void boundary_nodes_run(const Callback& callback) const {
        for(size_t b = 0; b < _mesh.boundary_groups_count(); ++b)
            for(size_t el = 0; el < _mesh.elements_count(b); ++el) {
                const auto& be = _mesh.element_1d(b, el);
                for(size_t i = 0; i < be->nodes_count(); ++i)
                    callback(b, el, i);
            }
    }

    void convert_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K, const std::vector<std::unordered_set<I>>& portrait) const {
        static constexpr auto accumulator = [](const size_t sum, const std::unordered_set<I>& row) { return sum + row.size(); };
        K.data().resize(std::accumulate(portrait.cbegin(), portrait.cend(), size_t{0}, accumulator));
        K.outerIndexPtr()[0] = 0;
        for(size_t row = 0; row < K.rows(); ++row) {
            K.outerIndexPtr()[row+1] = K.outerIndexPtr()[row];
            if (row < portrait.size())
                K.outerIndexPtr()[row+1] += portrait[row].size();
        }

#pragma omp parallel for default(none) shared(K, portrait)
        for(size_t row = 0; row < portrait.size(); ++row) {
            I inner_index = K.outerIndexPtr()[row];
            for(const I col : portrait[row]) {
                K.valuePtr()[inner_index] = 0;
                K.innerIndexPtr()[inner_index++] = col;
            }
            std::sort(&K.innerIndexPtr()[K.outerIndexPtr()[row]], &K.innerIndexPtr()[K.outerIndexPtr()[row+1]]);
        }
    }

    void approx_quad_nodes_on_bound(std::vector<std::array<T, 2>>& quad_nodes, const size_t b, const size_t el) const {
        const auto& be = _mesh.element_1d(b, el);
        quad_nodes.clear();
        quad_nodes.resize(be->qnodes_count(), {});
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            for(size_t i = 0; i < be->nodes_count(); ++i)
                for(size_t comp = 0; comp < 2; ++comp)
                    quad_nodes[q][comp] += _mesh.node(_mesh.node_number(b, el, i))[comp] * be->qN(i, q);
    }

    void approx_jacobi_matrices_on_bound(std::vector<std::array<T, 2>>& jacobi_matrices, const size_t b, const size_t el) const {
        const auto& be = _mesh.element_1d(b, el);
        jacobi_matrices.clear();
        jacobi_matrices.resize(be->qnodes_count(), {});
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            for(size_t i = 0; i < be->nodes_count(); ++i)
                for(size_t comp = 0; comp < 2; ++comp)
                    jacobi_matrices[q][comp] += _mesh.node(_mesh.node_number(b, el, i))[comp] * be->qNxi(i, q);
    }

    // Function - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Finite_Element_2D_Ptr, class Function>
    T integrate_function(const Finite_Element_2D_Ptr& e, const size_t i, size_t quad_shift, const Function& func) const {
        T integral = 0;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) * e->qN(i, q) * func(_quad_coords[quad_shift]) * jacobian(_jacobi_matrices[quad_shift]);
        return integral;
    }

    // Boundary_Gradient - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Boundary_Gradient>
    T integrate_boundary_gradient(const Finite_Element_1D_Ptr& be, const size_t i,
                                  const std::vector<std::array<T, 2>>& quad_nodes, 
                                  const std::vector<std::array<T, 2>>& jacobi_matrices, 
                                  const Boundary_Gradient& boundary_gradient) const {
        T integral = 0;
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            integral += be->weight(q) * be->qN(i, q) * boundary_gradient(quad_nodes[q]) * jacobian(jacobi_matrices[q]);
        return integral;
    }

    template<class B, size_t N>
    void integrate_boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const std::vector<boundary_condition<T, B, N>>& bounds_cond) const {
        std::vector<std::array<T, 2>> quad_nodes, jacobi_matrices;
        for(size_t b = 0; b < bounds_cond.size(); ++b)
            if (bounds_cond[b].contains_condition_second_kind())
                for(size_t el = 0; el < mesh().elements_count(b); ++el) {
                    approx_quad_nodes_on_bound(quad_nodes, b, el);
                    approx_jacobi_matrices_on_bound(jacobi_matrices, b, el);
                    const auto& be = mesh().element_1d(b, el);
                    for(size_t i = 0; i < be->nodes_count(); ++i)
                        for(size_t comp = 0; comp < bounds_cond[b].degrees_of_freedom(); ++comp)
                            if(bounds_cond[b].type(comp) == B(boundary_type::SECOND_KIND))
                                f[bounds_cond[b].degrees_of_freedom()*mesh().node_number(b, el, i) + comp] +=
                                    integrate_boundary_gradient(be, i, quad_nodes, jacobi_matrices, bounds_cond[b].func(comp));
                }
    }

    template<class B, size_t N>
    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const std::vector<boundary_condition<T, B, N>>& bounds_cond,
                                       const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound) const {
        Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(f.size());

        boundary_nodes_run(
            [this, &bounds_cond, &x](const size_t b, const size_t el, const size_t i) {
                for(size_t comp = 0; comp < bounds_cond[b].degrees_of_freedom(); ++comp)
                    if (bounds_cond[b].type(comp) == B(boundary_type::FIRST_KIND)) {
                        const I node = bounds_cond[b].degrees_of_freedom() * mesh().node_number(b, el, i) + comp;
                        if (x[node] == 0)
                            x[node] = bounds_cond[b].func(comp)(mesh().node(node));
                    }
            });

        f -= K_bound * x;

        boundary_nodes_run(
            [this, &bounds_cond, &x, &f](const size_t b, const size_t el, const size_t i) {
                for(size_t comp = 0; comp < bounds_cond[b].degrees_of_freedom(); ++comp)
                    if (bounds_cond[b].type(comp) == B(boundary_type::FIRST_KIND)) {
                        const I node = bounds_cond[b].degrees_of_freedom() * mesh().node_number(b, el, i) + comp;
                        f[node] = x[node];
                    }
            });
    }

    T jacobian(const size_t quad_shift) const noexcept {
        return jacobian(_jacobi_matrices[quad_shift]);
    }

    static T jacobian(const std::array<T, 4>& jacobi_matrix) noexcept {
        return std::abs(jacobi_matrix[0] * jacobi_matrix[3] - jacobi_matrix[1] * jacobi_matrix[2]);
    }

    static T jacobian(const std::array<T, 2>& jacobi_matrix) noexcept {
        return sqrt(jacobi_matrix[0] * jacobi_matrix[0] + jacobi_matrix[1] * jacobi_matrix[1]);
    }

    template<bool Component>
    T dNd(const Finite_Element_2D_Ptr& e, const size_t i, const size_t q, const size_t quad_shift) const {
        return dNd<Component>(e, i, q, _jacobi_matrices[quad_shift]); 
    }

    template<bool Component>
    static T dNd(const Finite_Element_2D_Ptr& e, const size_t i, const size_t q, const std::array<T, 4>& jacobi_matrix) {
        if constexpr (Component == component::X)
            return  e->qNxi(i, q) * jacobi_matrix[3] - e->qNeta(i, q) * jacobi_matrix[2];
        if constexpr (Component == component::Y)
            return -e->qNxi(i, q) * jacobi_matrix[1] + e->qNeta(i, q) * jacobi_matrix[0];
    }

public:
    void find_neighbors(const T r) {
        std::unordered_set<I> elements;
        for(size_t node = _first_node; node < _last_node; ++node)
            for(const I element : _nodes_elements_map[node])
                elements.insert(element);

        const std::vector<std::array<T, 2>> centres = approx_centres_of_elements(_mesh);
        _elements_neighbors = find_elements_neighbors(_mesh, centres, elements, r);

        double sum = 0;
        for(size_t i = 0; i < _elements_neighbors.size(); ++i)
            sum += _elements_neighbors[i].size();
        sum /= _elements_neighbors.size();
        std::cout << "Average neighbours count: " << sum << std::endl;
    }

    void set_neighbors(std::vector<std::vector<I>>&& neighbors) {
        _elements_neighbors = std::move(neighbors);
    }
};

}

#endif