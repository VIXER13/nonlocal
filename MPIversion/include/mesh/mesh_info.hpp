#ifndef MESH_INFO_HPP
#define MESH_INFO_HPP

#include "mesh_2d.hpp"
#include "utils.hpp"
#include <unordered_set>
#include <unordered_map>
#include <mpi.h>

namespace mesh {

enum class balancing_t : uint8_t { NO, MEMORY, SPEED };

template<class T, class I>
class mesh_info final {
    std::shared_ptr<mesh_2d<T, I>>              _mesh;
    std::vector<I>                              _quad_shifts;               // Квадратурные сдвиги.
    std::vector<std::array<T, 2>>               _quad_coords;               // Координаты квадратурных узлов сетки.
    std::vector<std::array<T, 4>>               _jacobi_matrices,           // Матрицы Якоби вычисленные в квадратурных узлах.
                                                _jacobi_matrices_nodes;     // Матрицы Якоби вычисленные в узлах сетки.
    std::vector<std::vector<I>>                 _quad_shifts_bound;         // Квадратурные сдвиги для граничных элементов.
    std::vector<std::vector<std::array<T, 2>>>  _quad_coords_bound,         // Координаты квадратурных узлов на границе.
                                                _jacobi_matrices_bound;     // Матрицы Якоби на границе.
    std::vector<std::vector<I>>                 _nodes_elements_map;        // Номера элементов, в которых присутствует узел
    std::vector<std::unordered_map<I, uint8_t>> _global_to_local_numbering; // Переход от глобальной нумерации к локальной каждого элемента
                                                                            // Считаем, что в элементе не более 255 узлов.
    std::vector<std::vector<I>>                 _elements_neighbors;        // Массив с номерами ближайших соседей
    int                                         _rank       = 0, _size      = 0;
    size_t                                      _first_node = 0, _last_node = 0;

    // Квадратурные сдвиги по элементам.
    static std::vector<I> quadrature_shifts_init(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        std::vector<I> quad_shifts(mesh->elements_count()+1);
        quad_shifts[0] = 0;
        for(size_t e = 0; e < mesh->elements_count(); ++e)
            quad_shifts[e+1] = quad_shifts[e] + mesh->element_2d(e)->qnodes_count();
        return std::move(quad_shifts);
    }

    static void check_shifts(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<I>& quad_shifts) {
        if(mesh->elements_count()+1 != quad_shifts.size())
            throw std::logic_error{"Quadrature shifts are incorrect size."};
    }

    // Аппроксимация глобальных координат всех квадратурных узлов сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    static std::vector<std::array<T, 2>> approx_all_quad_nodes(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<I>& quad_shifts) {
        check_shifts(mesh, quad_shifts);
        std::vector<std::array<T, 2>> quad_coords(quad_shifts.back(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, quad_coords)
        for(size_t e = 0; e < mesh->elements_count(); ++e) {
            const auto& el = mesh->element_2d(e);
            for(size_t q = 0; q < el->qnodes_count(); ++q)
                for(size_t i = 0; i < el->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 2; ++comp)
                        quad_coords[quad_shifts[e]+q][comp] += mesh->node(mesh->node_number(e, i))[comp] * el->qN(i, q);
        }
        return std::move(quad_coords);
    }

    // Аппроксимация матриц Якоби во всех квадратурных узлах сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    static std::vector<std::array<T, 4>> approx_all_jacobi_matrices(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<I>& quad_shifts) {
        check_shifts(mesh, quad_shifts);
        std::vector<std::array<T, 4>> jacobi_matrices(quad_shifts.back(), std::array<T, 4>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, jacobi_matrices)
        for(size_t e = 0; e < mesh->elements_count(); ++e) {
            const auto& el = mesh->element_2d(e);
            for(size_t q = 0; q < el->qnodes_count(); ++q)
                for(size_t i = 0; i < el->nodes_count(); ++i) {
                    jacobi_matrices[quad_shifts[e]+q][0] += mesh->node(mesh->node_number(e, i))[0] * el->qNxi (i, q);
                    jacobi_matrices[quad_shifts[e]+q][1] += mesh->node(mesh->node_number(e, i))[0] * el->qNeta(i, q);
                    jacobi_matrices[quad_shifts[e]+q][2] += mesh->node(mesh->node_number(e, i))[1] * el->qNxi (i, q);
                    jacobi_matrices[quad_shifts[e]+q][3] += mesh->node(mesh->node_number(e, i))[1] * el->qNeta(i, q);
                }
        }
        return std::move(jacobi_matrices);
    }

    std::vector<std::array<T, 4>> approx_jacobi_matrices_nodes() {
        std::vector<std::array<T, 4>> jacobi_matrices(mesh().nodes_count(), std::array<T, 4>{});
#pragma omp parallel for default(none) shared(jacobi_matrices)
        for(size_t node = 0; node < mesh().nodes_count(); ++node) {
            for(const I e : nodes_elements_map(node)) {
                const auto& el = mesh().element_2d(e);
                for(size_t i = 0; i < el->nodes_count(); ++i) {
                    for(size_t j = 0; j < el->nodes_count(); ++j) {
                        const std::array<T, 2>& mesh_node = mesh().node(mesh().node_number(e, j));
                        const T Nxi  = el->Nxi (j, el->node(i)),
                                Neta = el->Neta(j, el->node(i));
                        jacobi_matrices[node][0] += mesh_node[0] * Nxi;
                        jacobi_matrices[node][1] += mesh_node[0] * Neta;
                        jacobi_matrices[node][2] += mesh_node[1] * Nxi;
                        jacobi_matrices[node][3] += mesh_node[1] * Neta;
                    }
                }
            }
            for(size_t i = 0; i < 4; ++i)
                jacobi_matrices[node][i] /= nodes_elements_map(node).size();
        }
        return std::move(jacobi_matrices);
    }

    // Квадратурные сдвиги для граничных элементов.
    static std::vector<std::vector<I>> quadrature_shifts_bound_init(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        std::vector<std::vector<I>> quad_shifts(mesh->boundary_groups_count());
        for(size_t b = 0; b < mesh->boundary_groups_count(); ++b) {
            quad_shifts[b].resize(mesh->elements_count(b) + 1);
            quad_shifts[b][0] = 0;
            for(size_t e = 0; e < mesh->elements_count(b); ++e)
                quad_shifts[b][e+1] = quad_shifts[b][e] + mesh->element_1d(b, e)->qnodes_count();
        }
        return std::move(quad_shifts);
    }

    static void check_shifts(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<std::vector<I>>& quad_shifts) {
        if(mesh->boundary_groups_count() != quad_shifts.size())
            throw std::logic_error{"The number of boundary groups does not match the number of shifts groups."};
        for(size_t b = 0; b < mesh->boundary_groups_count(); ++b)
            if(mesh->elements_count(b)+1 != quad_shifts[b].size())
                throw std::logic_error{"Quadrature shifts on bound are incorrect size."};
    }

    // Координаты квадратурных узлов на границах.
    static std::vector<std::vector<std::array<T, 2>>> approx_all_quad_nodes_bound(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<std::vector<I>>& quad_shifts) {
        check_shifts(mesh, quad_shifts);
        std::vector<std::vector<std::array<T, 2>>> quad_coords(mesh->boundary_groups_count());
        for(size_t b = 0; b < quad_coords.size(); ++b) {
            quad_coords[b].resize(quad_shifts[b].back(), std::array<T, 2>{});
            for(size_t e = 0; e < mesh->elements_count(b); ++e) {
                const auto& el = mesh->element_1d(b, e);
                for(size_t q = 0; q < el->qnodes_count(); ++q)
                    for(size_t i = 0; i < el->nodes_count(); ++i)
                        for(size_t comp = 0; comp < 2; ++comp)
                            quad_coords[b][quad_shifts[b][e]+q][comp] += mesh->node(mesh->node_number(b, e, i))[comp] * el->qN(i, q);
            }
        }
        return std::move(quad_coords);
    }

    // Матрицы якоби на границах.
    static std::vector<std::vector<std::array<T, 2>>> approx_all_jacobi_matrices_bound(const std::shared_ptr<mesh_2d<T, I>>& mesh, const std::vector<std::vector<I>>& quad_shifts) {
        check_shifts(mesh, quad_shifts);
        std::vector<std::vector<std::array<T, 2>>> jacobi_matrices(mesh->boundary_groups_count());
        for(size_t b = 0; b < jacobi_matrices.size(); ++b) {
            jacobi_matrices[b].resize(quad_shifts[b].back(), std::array<T, 2>{});
            for(size_t e = 0; e < mesh->elements_count(b); ++e) {
                const auto& el = mesh->element_1d(b, e);
                for(size_t q = 0; q < el->qnodes_count(); ++q)
                    for(size_t i = 0; i < el->nodes_count(); ++i)
                        for(size_t comp = 0; comp < 2; ++comp)
                            jacobi_matrices[b][quad_shifts[b][e]+q][comp] += mesh->node(mesh->node_number(b, e, i))[comp] * el->qNxi(i, q);
            }
        }
        return std::move(jacobi_matrices);
    }

    static std::vector<std::vector<I>> node_elements_map_init(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        std::vector<std::vector<I>> nodes_elements_map(mesh->nodes_count());
        for(size_t e = 0; e < mesh->elements_count(); ++e) {
            const auto& el = mesh->element_2d(e);
            for(size_t node = 0; node < el->nodes_count(); ++node)
                nodes_elements_map[mesh->node_number(e, node)].push_back(e);
        }
        return std::move(nodes_elements_map);
    }

    static std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering_init(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh->elements_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
        for(size_t e = 0; e < mesh->elements_count(); ++e) {
            const auto& el = mesh->element_2d(e);
            for(size_t node = 0; node < el->nodes_count(); ++node)
                global_to_local_numbering[e][mesh->node_number(e, node)] = node;
        }
        return std::move(global_to_local_numbering);
    }

    static std::vector<std::array<T, 2>> approx_centres_of_elements(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        std::vector<std::array<T, 2>> centres(mesh->elements_count(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, centres)
        for(size_t e = 0; e < centres.size(); ++e) {
            const auto& el = mesh->element_2d(e);
            const T x0 = is_trinagle(mesh->element_2d_type(e)) ? T{1}/T{3} : T{0};
            for(size_t node = 0; node < el->nodes_count(); ++node) {
                using namespace utils;
                centres[e] += mesh->node(mesh->node_number(e, node)) * el->N(node, {x0, x0});
            }
        }
        return std::move(centres);
    }

    static void find_elements_neighbors(std::vector<std::vector<I>>& elements_neighbors,
                                        const std::shared_ptr<mesh_2d<T, I>>& mesh,
                                        const std::vector<std::array<T, 2>>& centres,
                                        const std::unordered_set<I>& elements, const T r) {
        elements_neighbors.resize(mesh->elements_count());
        for(const I elL : elements)
            if (elements_neighbors[elL].empty()) {
                elements_neighbors[elL].reserve(mesh->elements_count());
                for(size_t elNL = 0; elNL < mesh->elements_count(); ++elNL)
                    if(utils::distance(centres[elL], centres[elNL]) < r)
                        elements_neighbors[elL].push_back(elNL);
                elements_neighbors[elL].shrink_to_fit();
            }
    }

public:
    explicit mesh_info(const std::shared_ptr<mesh_2d<T, I>>& mesh) { set_mesh(mesh); }

    void set_mesh(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
        if (mesh == nullptr)
            throw std::invalid_argument{"mesh can't nullptr"};
        _mesh = mesh;
        _quad_shifts = quadrature_shifts_init(_mesh);
        _quad_coords = approx_all_quad_nodes(_mesh, _quad_shifts);
        _jacobi_matrices = approx_all_jacobi_matrices(_mesh, _quad_shifts);
        _quad_shifts_bound = quadrature_shifts_bound_init(_mesh);
        _quad_coords_bound = approx_all_quad_nodes_bound(_mesh, _quad_shifts_bound);
        _jacobi_matrices_bound = approx_all_jacobi_matrices_bound(_mesh, _quad_shifts_bound);
        _nodes_elements_map = node_elements_map_init(_mesh);
        _global_to_local_numbering = global_to_local_numbering_init(_mesh);
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        _first_node = _mesh->nodes_count() / _size *  _rank;
        _last_node  = _mesh->nodes_count() / _size * (_rank+1) + (_rank == _size-1) * _mesh->nodes_count() % _size;

        _jacobi_matrices_nodes = approx_jacobi_matrices_nodes();
    }

    const mesh_2d<T, I>&                  mesh                     ()                     const { return *_mesh; }
    int                                   rank                     ()                     const { return _rank; }
    int                                   size                     ()                     const { return _size; }
    size_t                                first_node               ()                     const { return _first_node; }
    size_t                                last_node                ()                     const { return _last_node; }
    I                                     quad_shift               (const size_t element) const { return _quad_shifts[element]; }
    const std::array<T, 2>&               quad_coord               (const size_t quad)    const { return _quad_coords[quad]; }
    const std::array<T, 4>&               jacobi_matrix            (const size_t quad)    const { return _jacobi_matrices[quad]; }
    const std::array<T, 4>&               jacobi_matrix_node       (const size_t node)    const { return _jacobi_matrices_nodes[node]; }
    const std::vector<I>&                 nodes_elements_map       (const size_t node)    const { return _nodes_elements_map[node]; }
    const std::unordered_map<I, uint8_t>& global_to_local_numbering(const size_t element) const { return _global_to_local_numbering[element]; }
    const std::vector<I>&                 neighbors                (const size_t element) const { return _elements_neighbors[element]; }
    I                                     quad_shift               (const size_t bound, const size_t element) const { return _quad_shifts_bound[bound][element]; }
    const std::array<T, 2>&               quad_coord               (const size_t bound, const size_t quad)    const { return _quad_coords_bound[bound][quad]; }
    const std::array<T, 2>&               jacobi_matrix            (const size_t bound, const size_t quad)    const { return _jacobi_matrices_bound[bound][quad]; }

    static T jacobian(const std::array<T, 4>& J) noexcept { return std::abs(J[0] * J[3] - J[1] * J[2]); }
           T jacobian(const size_t quad)         const    { return jacobian(jacobi_matrix(quad)); }

    void find_neighbours(const T r, const balancing_t balancing) {
        _elements_neighbors.clear();
        const std::vector<std::array<T, 2>> centres = approx_centres_of_elements(_mesh);
        std::unordered_set<I> elements;
        for(size_t node = first_node(); node < last_node(); ++node)
            for(const I e : nodes_elements_map(node))
                elements.insert(e);
        find_elements_neighbors(_elements_neighbors, _mesh, centres, elements, r);

        if (size() > 1 && balancing != balancing_t::NO) {
            std::vector<int> data_count_per_nodes(mesh().nodes_count(), 0);
            switch(balancing) {
                case balancing_t::MEMORY:
                    elements.clear();
#pragma omp parallel for default(none) shared(data_count_per_nodes) private(elements)
                    for(size_t node = first_node(); node < last_node(); ++node)
                        for(const I e : nodes_elements_map(node)) {
                            const auto [it, inserted] = elements.insert(e);
                            if (inserted)
                                data_count_per_nodes[node] += neighbors(e).size();
                        }
                break;

                case balancing_t::SPEED:
#pragma omp parallel for default(none) shared(data_count_per_nodes)
                    for(size_t node = first_node(); node < last_node(); ++node)
                        for(const I eL : nodes_elements_map(node))
                            for(const I eNL : neighbors(eL))
                                for(size_t i = 0; i < mesh().element_2d(eNL)->nodes_count(); ++i)
                                    if (node <= mesh().node_number(eNL, i))
                                        ++data_count_per_nodes[node];
                break;
            }

            const std::vector<int> sendcounts(size(), last_node() - first_node()), sdispls(size(), first_node());
            std::vector<int> recvcounts(size()), rdispls(size());
            std::vector<std::array<size_t, 2>> first_last(size());
            for(int i = 0; i < size(); ++i) {
                first_last[i].front() = mesh().nodes_count() / size() *  i;
                first_last[i].back()  = mesh().nodes_count() / size() * (i+1) + (i == size()-1) * mesh().nodes_count() % size();
                recvcounts[i] = first_last[i].back() - first_last[i].front();
                rdispls[i]    = first_last[i].front();
            }
            MPI_Alltoallv(data_count_per_nodes.data(), sendcounts.data(), sdispls.data(), MPI_INT,
                          data_count_per_nodes.data(), recvcounts.data(), rdispls.data(), MPI_INT, MPI_COMM_WORLD);

            size_t sum = 0, curr_rank = 0;
            const size_t mean = std::accumulate(data_count_per_nodes.cbegin(), data_count_per_nodes.cend(), size_t{0}) / size();
            for(size_t node = 0; node < data_count_per_nodes.size(); ++node) {
                sum += data_count_per_nodes[node];
                if (sum > mean) {
                    first_last[curr_rank].back() = node;
                    first_last[curr_rank+1].front() = node;
                    ++curr_rank;
                    sum = 0;
                }
            }
            _first_node = first_last[rank()].front();
            _last_node  = first_last[rank()].back();

            elements.clear();
            for(size_t node = first_node(); node < last_node(); ++node)
                for(const I e : nodes_elements_map(node))
                    elements.insert(e);
            find_elements_neighbors(_elements_neighbors, _mesh, centres, elements, r);
            for(size_t e = 0; e < _elements_neighbors.size(); ++e)
                if (elements.find(e) == elements.cend()) {
                    _elements_neighbors[e].clear();
                    _elements_neighbors[e].shrink_to_fit();
                }
        }
    }

    T integrate_solution(const std::vector<T>& sol) const {
        if(mesh().nodes_count() != size_t(sol.size()))
            throw std::logic_error{"mesh.nodes_count() != T.size()"};

        T integral = 0;
        for(size_t el = 0; el < mesh().elements_count(); ++el) {
            const auto& e = mesh().element_2d(el);
            for(size_t i = 0; i < e->nodes_count(); ++i)
                for(size_t q = 0, shift = quad_shift(el); q < e->qnodes_count(); ++q, ++shift)
                    integral += e->weight(q) * e->qN(i, q) * sol[mesh().node_number(el, i)] * jacobian(shift);
        }
        return integral;
    }

    std::array<std::vector<T>, 2> calc_gradient(const std::vector<T>& sol) const {
        if(mesh().nodes_count() != sol.size())
            throw std::logic_error{"mesh.nodes_count() != sol.size()"};

        std::vector<T> x(sol.size(), 0), y(sol.size(), 0);
#pragma omp parallel for default(none) shared(x, y, sol)
        for(size_t node = 0; node < mesh().nodes_count(); ++node) {
            T dx = 0, dy = 0;
            const T jac = jacobian(jacobi_matrix_node(node));
            for(const I e : nodes_elements_map(node)) {
                const auto& el = mesh().element_2d(e);
                for(size_t i = 0; i < el->nodes_count(); ++i)
                    for(size_t j = 0; j < el->nodes_count(); ++j) {
                        const T Nxi  = el->Nxi (j, el->node(i)) * sol[mesh().node_number(e, j)] / jac,
                                Neta = el->Neta(j, el->node(i)) * sol[mesh().node_number(e, j)] / jac;
                        dx +=  jacobi_matrix_node(node)[3] * Nxi - jacobi_matrix_node(node)[2] * Neta;
                        dy += -jacobi_matrix_node(node)[1] * Nxi + jacobi_matrix_node(node)[0] * Neta;
                    }
            }
            x[node] += dx / nodes_elements_map(node).size();
            y[node] += dy / nodes_elements_map(node).size();
        }
        return {std::move(x), std::move(y)};
    }

    std::vector<T> approx_in_quad(const std::vector<T>& x) {
        std::vector<T> x_in_quad(quad_shift(mesh().elements_count()), 0);
#pragma omp parallel for default(none) shared(x, x_in_quad)
        for(size_t el = 0; el < mesh().elements_count(); ++el) {
            const auto& e = mesh().element_2d(el);
            for(size_t q = 0, shift = quad_shift(el); q < e->qnodes_count(); ++q, ++shift)
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    x_in_quad[shift] += x[mesh().node_number(el, i)] * e->qN(i, q);
        }
        return std::move(x_in_quad);
    }
};

}

#endif