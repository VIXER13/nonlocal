#ifndef MESH_PROXY_HPP
#define MESH_PROXY_HPP

#include "mesh_container_2d.hpp"
#include "utils.hpp"
#include "MPI_utils.hpp"
#include <unordered_set>
#include <unordered_map>
#include <numeric>

namespace nonlocal::mesh {

enum class balancing_t : uint8_t { NO, MEMORY, SPEED };

template<class T, class I>
class mesh_proxy final {
    std::shared_ptr<mesh_2d<T, I>>              _mesh;                      // Сетка из файла.
    std::vector<std::vector<I>>                 _nodes_elements_map;        // Номера элементов, в которых присутствует узел.
    std::vector<std::unordered_map<I, uint8_t>> _global_to_local_numbering; // Переход от глобальной нумерации к локальной каждого элемента.
                                                                            // Считаем, что в элементе не более 255 узлов.

    std::vector<I>                              _quad_shifts;               // Квадратурные сдвиги.
    std::vector<std::array<T, 2>>               _quad_coords;               // Координаты квадратурных узлов сетки.
    std::vector<std::array<T, 4>>               _jacobi_matrices;           // Матрицы Якоби вычисленные в квадратурных узлах.

    std::vector<I>                              _quad_node_shift;           // Квадратурные сдвиги помноженные на количество узлов.
    std::vector<std::array<T, 2>>               _dNdX;                      // Массив с вычисленными глобальными производными

    std::unordered_map<std::string_view, std::vector<I>>                _quad_shifts_bound;     // Квадратурные сдвиги для граничных элементов.
    std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> _quad_coords_bound,     // Координаты квадратурных узлов на границе.
                                                                        _jacobi_matrices_bound; // Матрицы Якоби на границе.
    std::vector<T>                              _elements_ares;
    std::vector<std::array<T, 4>>               _jacobi_matrices_nodes;     // Матрицы Якоби вычисленные в узлах сетки.

    MPI_utils::MPI_ranges                       _ranges;                    // С какого и по какой узлы распределена обработка данных между процессами

    std::vector<std::vector<I>>                 _elements_neighbors;        // Массив с номерами ближайших соседей.

    static std::vector<std::vector<I>>                 node_elements_map_init        (const mesh_2d<T, I>& mesh);
    static std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering_init(const mesh_2d<T, I>& mesh);

    static std::vector<I>                quadrature_shifts_init    (const mesh_2d<T, I>& mesh);
    static void                          check_shifts              (const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts);
    static std::vector<std::array<T, 2>> approx_all_quad_nodes     (const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts);
    static std::vector<std::array<T, 4>> approx_all_jacobi_matrices(const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts);

    static std::vector<I> quadrature_node_shits_init(const mesh_2d<T, I>& mesh);
    static std::vector<std::array<T, 2>> dNdX_init  (const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts,
                                                     const std::vector<std::array<T, 4>>& jacobi_matrices,
                                                     const std::vector<I>& quad_element_shift);

    static std::unordered_map<std::string_view, std::vector<I>>                quadrature_shifts_bound_init    (const mesh_2d<T, I>& mesh);
    static void                                                                check_shifts                    (const mesh_2d<T, I>& mesh, const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts);
    static std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> approx_all_quad_nodes_bound     (const mesh_2d<T, I>& mesh, const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts);
    static std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> approx_all_jacobi_matrices_bound(const mesh_2d<T, I>& mesh, const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts);

    static std::vector<T> approx_elements_areas(const mesh_2d<T, I>& mesh,
                                                const std::vector<I>& quad_shifts,
                                                const std::vector<std::array<T, 4>>& jacobi_matrices);
    static std::vector<std::array<T, 4>> approx_jacobi_matrices_nodes(const mesh_2d<T, I>& mesh,
                                                                      const std::vector<T>& elements_ares,
                                                                      const std::vector<std::vector<I>>& nodes_elements_map,
                                                                      const std::vector<std::unordered_map<I, uint8_t>>& global_to_local_numbering);

    void first_last_node_init();

    static std::vector<std::array<T, 2>> approx_centres_of_elements(const mesh_2d<T, I>& mesh);
    static void find_elements_neighbors(std::vector<std::vector<I>>& elements_neighbors,
                                        const mesh_2d<T, I>& mesh,
                                        const std::vector<std::array<T, 2>>& centres,
                                        const std::unordered_set<I>& elements, const T r);

    I quad_shift(const std::string_view bound, const size_t element) const;
    I quad_node_shift(const size_t element, const size_t i) const;

public:
    using quad_coord_iterator          = typename std::vector<std::array<T, 2>>::const_iterator;
    using jacobi_matrix_iterator       = typename std::vector<std::array<T, 4>>::const_iterator;
    using dNdX_iterator                = typename std::vector<std::array<T, 2>>::const_iterator;
    using quad_coord_bound_iterator    = typename std::vector<std::array<T, 2>>::const_iterator;
    using jacobi_matrix_bound_iterator = typename std::vector<std::array<T, 2>>::const_iterator;

    explicit mesh_proxy(const std::shared_ptr<mesh_2d<T, I>>& mesh);

    void set_mesh(const std::shared_ptr<mesh_2d<T, I>>& mesh);

    const mesh_2d<T, I>&                  mesh                     ()                                         const;
    const std::vector<I>&                 nodes_elements_map       (const size_t node)                        const;
    const std::unordered_map<I, uint8_t>& global_to_local_numbering(const size_t element)                     const;
    const std::array<T, 4>&               jacobi_matrix_node       (const size_t node)                        const;
    T                                     element_area             (const size_t element)                     const;
    I                                     quad_shift               (const size_t element)                     const;
    quad_coord_iterator                   quad_coord               (const size_t element)                     const;
    jacobi_matrix_iterator                jacobi_matrix            (const size_t element)                     const;
    dNdX_iterator                         dNdX                     (const size_t element, const size_t i)     const;
    quad_coord_bound_iterator             quad_coord               (const std::string_view bound, const size_t element) const;
    jacobi_matrix_bound_iterator          jacobi_matrix            (const std::string_view bound, const size_t element) const;
    MPI_utils::MPI_ranges                 ranges                   ()                                         const;
    size_t                                first_node               ()                                         const;
    size_t                                last_node                ()                                         const;
    const std::vector<I>&                 neighbors                (const size_t element)                     const;

    static T jacobian(const std::array<T, 2>& J) noexcept;
    static T jacobian(const std::array<T, 4>& J) noexcept;

    void find_neighbours(const T r, const balancing_t balancing);

    template<class Vector>
    T integrate_solution(const Vector& sol) const;
    template<class Vector>
    std::array<std::vector<T>, 2> gradient(const Vector& sol) const;
    std::vector<T> approx_in_quad(const std::vector<T>& x);
};

template<class T, class I>
mesh_proxy<T, I>::mesh_proxy(const std::shared_ptr<mesh_2d<T, I>>& mesh) { set_mesh(mesh); }

template<class T, class I>
void mesh_proxy<T, I>::set_mesh(const std::shared_ptr<mesh_2d<T, I>>& mesh) {
    if (mesh == nullptr)
        throw std::invalid_argument{"mesh can't nullptr"};

    _mesh = mesh;
    _nodes_elements_map = node_elements_map_init(*_mesh);
    _global_to_local_numbering = global_to_local_numbering_init(*_mesh);

    _quad_shifts = quadrature_shifts_init(*_mesh);
    _quad_coords = approx_all_quad_nodes(*_mesh, _quad_shifts);
    _jacobi_matrices = approx_all_jacobi_matrices(*_mesh, _quad_shifts);

    _quad_node_shift = quadrature_node_shits_init(*_mesh);
    _dNdX = dNdX_init(*_mesh, _quad_shifts, _jacobi_matrices, _quad_node_shift);

    _quad_shifts_bound = quadrature_shifts_bound_init(*_mesh);
    _quad_coords_bound = approx_all_quad_nodes_bound(*_mesh, _quad_shifts_bound);
    _jacobi_matrices_bound = approx_all_jacobi_matrices_bound(*_mesh, _quad_shifts_bound);

    _elements_ares = approx_elements_areas(*_mesh, _quad_shifts, _jacobi_matrices);
    _jacobi_matrices_nodes = approx_jacobi_matrices_nodes(*_mesh, _elements_ares, _nodes_elements_map, _global_to_local_numbering);

    first_last_node_init();
}

template<class T, class I>
const mesh_2d<T, I>& mesh_proxy<T, I>::mesh() const { return *_mesh; }

template<class T, class I>
const std::vector<I>& mesh_proxy<T, I>::nodes_elements_map(const size_t node) const { return _nodes_elements_map[node]; }

template<class T, class I>
const std::unordered_map<I, uint8_t>& mesh_proxy<T, I>::global_to_local_numbering(const size_t element) const { return _global_to_local_numbering[element]; }

template<class T, class I>
const std::array<T, 4>& mesh_proxy<T, I>::jacobi_matrix_node(const size_t node) const { return _jacobi_matrices_nodes[node]; }

template<class T, class I>
T mesh_proxy<T, I>::element_area(const size_t element) const { return _elements_ares[element]; }

template<class T, class I>
I mesh_proxy<T, I>::quad_shift(const size_t element) const { return _quad_shifts[element]; }

template<class T, class I>
typename mesh_proxy<T, I>::quad_coord_iterator mesh_proxy<T, I>::quad_coord(const size_t element) const {
    return std::next(_quad_coords.cbegin(), quad_shift(element));
}

template<class T, class I>
typename mesh_proxy<T, I>::jacobi_matrix_iterator mesh_proxy<T, I>::jacobi_matrix(const size_t element) const {
    return std::next(_jacobi_matrices.cbegin(), quad_shift(element));
}

template<class T, class I>
I mesh_proxy<T, I>::quad_shift(const std::string_view bound, const size_t element) const { return _quad_shifts_bound.at(bound)[element]; }

template<class T, class I>
typename mesh_proxy<T, I>::quad_coord_bound_iterator mesh_proxy<T, I>::quad_coord(const std::string_view bound, const size_t element) const {
    return std::next(_quad_coords_bound.at(bound).cbegin(), quad_shift(bound, element));
}

template<class T, class I>
typename mesh_proxy<T, I>::jacobi_matrix_bound_iterator mesh_proxy<T, I>::jacobi_matrix(const std::string_view bound, const size_t element) const {
    return std::next(_jacobi_matrices_bound.at(bound).cbegin(), quad_shift(bound, element));
}

template<class T, class I>
I mesh_proxy<T, I>::quad_node_shift(const size_t element, const size_t i) const {
    return _quad_node_shift[element] + i * mesh().element_2d(element)->qnodes_count();
}

template<class T, class I>
typename mesh_proxy<T, I>::dNdX_iterator mesh_proxy<T, I>::dNdX(const size_t element, const size_t i) const {
    return std::next(_dNdX.cbegin(), quad_node_shift(element, i));
}

template<class T, class I>
MPI_utils::MPI_ranges mesh_proxy<T, I>::ranges() const { return _ranges; }

template<class T, class I>
size_t mesh_proxy<T, I>::first_node() const { return _ranges.range().front(); }

template<class T, class I>
size_t mesh_proxy<T, I>::last_node() const { return _ranges.range().back(); }

template<class T, class I>
const std::vector<I>& mesh_proxy<T, I>::neighbors(const size_t element) const { return _elements_neighbors[element]; }

template<class T, class I>
T mesh_proxy<T, I>::jacobian(const std::array<T, 2>& J) noexcept { return std::sqrt(J[0] * J[0] + J[1] * J[1]); }

template<class T, class I>
T mesh_proxy<T, I>::jacobian(const std::array<T, 4>& J) noexcept { return std::abs(J[0] * J[3] - J[1] * J[2]); }

template<class T, class I>
std::vector<std::vector<I>> mesh_proxy<T, I>::node_elements_map_init(const mesh_2d<T, I>& mesh) {
    std::vector<std::vector<I>> nodes_elements_map(mesh.nodes_count());
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t node = 0; node < el->nodes_count(); ++node)
            nodes_elements_map[mesh.node_number(e, node)].push_back(e);
    }
    return nodes_elements_map;
}

template<class T, class I>
std::vector<std::unordered_map<I, uint8_t>> mesh_proxy<T, I>::global_to_local_numbering_init(const mesh_2d<T, I>& mesh) {
    std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh.elements_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t node = 0; node < el->nodes_count(); ++node)
            global_to_local_numbering[e][mesh.node_number(e, node)] = node;
    }
    return global_to_local_numbering;
}

template<class T, class I>
std::vector<T> mesh_proxy<T, I>::approx_elements_areas(const mesh_2d<T, I>& mesh,
                                                       const std::vector<I>& quad_shifts,
                                                       const std::vector<std::array<T, 4>>& jacobi_matrices
                                                       ) {
    std::vector<T> elements_areas(mesh.elements_count(), 0);
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q) {
            const T weight = el->weight(q) * jacobian(jacobi_matrices[quad_shifts[e]+q]);
            for(size_t i = 0; i < el->nodes_count(); ++i)
                elements_areas[e] += weight * el->qN(i, q);
        }
    }
    return elements_areas;
}

template<class T, class I>
std::vector<std::array<T, 4>> mesh_proxy<T, I>::approx_jacobi_matrices_nodes(const mesh_2d<T, I>& mesh,
                                                                             const std::vector<T>& elements_ares,
                                                                             const std::vector<std::vector<I>>& nodes_elements_map,
                                                                             const std::vector<std::unordered_map<I, uint8_t>>& global_to_local_numbering) {
    std::vector<std::array<T, 4>> jacobi_matrices(mesh.nodes_count(), std::array<T, 4>{});
#pragma omp parallel for default(none) shared(jacobi_matrices, mesh, elements_ares, nodes_elements_map, global_to_local_numbering)
    for(size_t node = 0; node < mesh.nodes_count(); ++node) {
        T node_area = T{0};
        for(const I e : nodes_elements_map[node]) {
            const auto& el = mesh.element_2d(e);
            const size_t i = global_to_local_numbering[e].find(node)->second;
            for(size_t j = 0; j < el->nodes_count(); ++j) {
                const std::array<T, 2>& mesh_node = mesh.node(mesh.node_number(e, j));
                const T Nxi  = el->Nxi (j, el->node(i)) * elements_ares[e], // TODO: optimize Nxi and Neta computations
                        Neta = el->Neta(j, el->node(i)) * elements_ares[e];
                jacobi_matrices[node][0] += mesh_node[0] * Nxi;
                jacobi_matrices[node][1] += mesh_node[0] * Neta;
                jacobi_matrices[node][2] += mesh_node[1] * Nxi;
                jacobi_matrices[node][3] += mesh_node[1] * Neta;
            }
            node_area += elements_ares[e];
        }
        using namespace metamath::function;
        jacobi_matrices[node] /= node_area;
    }
    return jacobi_matrices;
}

template<class T, class I>
std::vector<I> mesh_proxy<T, I>::quadrature_shifts_init(const mesh_2d<T, I>& mesh) {
    std::vector<I> quad_shifts(mesh.elements_count()+1);
    quad_shifts[0] = 0;
    for(size_t e = 0; e < mesh.elements_count(); ++e)
        quad_shifts[e+1] = quad_shifts[e] + mesh.element_2d(e)->qnodes_count();
    return quad_shifts;
}

template<class T, class I>
void mesh_proxy<T, I>::check_shifts(const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts) {
    if(mesh.elements_count()+1 != quad_shifts.size())
        throw std::logic_error{"Quadrature shifts are incorrect size."};
}

template<class T, class I>
std::vector<std::array<T, 2>> mesh_proxy<T, I>::approx_all_quad_nodes(const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts) {
    check_shifts(mesh, quad_shifts);
    std::vector<std::array<T, 2>> quad_coords(quad_shifts.back(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, quad_coords)
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q)
            for(size_t i = 0; i < el->nodes_count(); ++i)
                for(size_t comp = 0; comp < 2; ++comp)
                    quad_coords[quad_shifts[e]+q][comp] += mesh.node(mesh.node_number(e, i))[comp] * el->qN(i, q);
    }
    return quad_coords;
}

template<class T, class I>
std::vector<std::array<T, 4>> mesh_proxy<T, I>::approx_all_jacobi_matrices(const mesh_2d<T, I>& mesh, const std::vector<I>& quad_shifts) {
    check_shifts(mesh, quad_shifts);
    std::vector<std::array<T, 4>> jacobi_matrices(quad_shifts.back(), std::array<T, 4>{});
#pragma omp parallel for default(none) shared(mesh, quad_shifts, jacobi_matrices)
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q)
            for(size_t i = 0; i < el->nodes_count(); ++i) {
                jacobi_matrices[quad_shifts[e]+q][0] += mesh.node(mesh.node_number(e, i))[0] * el->qNxi (i, q);
                jacobi_matrices[quad_shifts[e]+q][1] += mesh.node(mesh.node_number(e, i))[0] * el->qNeta(i, q);
                jacobi_matrices[quad_shifts[e]+q][2] += mesh.node(mesh.node_number(e, i))[1] * el->qNxi (i, q);
                jacobi_matrices[quad_shifts[e]+q][3] += mesh.node(mesh.node_number(e, i))[1] * el->qNeta(i, q);
            }
    }
    return jacobi_matrices;
}

template<class T, class I>
std::vector<I> mesh_proxy<T, I>::quadrature_node_shits_init(const mesh_2d<T, I>& mesh) {
    std::vector<I> quad_element_shift(mesh.elements_count()+1);
    quad_element_shift[0] = 0;
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        quad_element_shift[e+1] = quad_element_shift[e] + el->nodes_count() * el->qnodes_count();
    }
    return quad_element_shift;
}

template<class T, class I>
std::vector<std::array<T, 2>> mesh_proxy<T, I>::dNdX_init(const mesh_2d<T, I>& mesh,
                                                          const std::vector<I>& quad_shifts,
                                                          const std::vector<std::array<T, 4>>& jacobi_matrices,
                                                          const std::vector<I>& quad_element_shift) {
    std::vector<std::array<T, 2>> dNdX(quad_element_shift.back());
    for(size_t e = 0; e < mesh.elements_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(size_t i = 0; i < el->nodes_count(); ++i)
            for(size_t q = 0; q < el->qnodes_count(); ++q) {
                const std::array<T, 4>& J = jacobi_matrices[quad_shifts[e] + q];
                dNdX[quad_element_shift[e] + i * el->qnodes_count() + q] = {
                        el->qNxi(i, q) * J[3] - el->qNeta(i, q) * J[2],
                        -el->qNxi(i, q) * J[1] + el->qNeta(i, q) * J[0]
                };
            }
    }
    return dNdX;
}

template<class T, class I>
std::unordered_map<std::string_view, std::vector<I>> mesh_proxy<T, I>::quadrature_shifts_bound_init(const mesh_2d<T, I>& mesh) {
    std::unordered_map<std::string_view, std::vector<I>> quad_shifts(mesh.boundary_groups_count());
    for(const std::string& bound : mesh.boundary_names()) {
        auto& shifts = quad_shifts[bound];
        shifts.resize(mesh.elements_count(bound) + 1);
        shifts[0] = 0;
        for(size_t e = 0; e < mesh.elements_count(bound); ++e)
            shifts[e+1] = shifts[e] + mesh.element_1d(bound, e)->qnodes_count();
    }
    return quad_shifts;
}

template<class T, class I>
void mesh_proxy<T, I>::check_shifts(const mesh_2d<T, I>& mesh, const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts) {
    if(mesh.boundary_groups_count() != quad_shifts.size())
        throw std::logic_error{"The number of boundary groups does not match the number of shifts groups."};
    for(const std::string& bound : mesh.boundary_names())
        if(mesh.elements_count(bound)+1 != quad_shifts.at(bound).size())
            throw std::logic_error{"Quadrature shifts on bound are incorrect size."};
}

template<class T, class I>
std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> mesh_proxy<T, I>::approx_all_quad_nodes_bound(const mesh_2d<T, I>& mesh,
                                                                                                                  const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts) {
    check_shifts(mesh, quad_shifts);
    std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> quad_coords(mesh.boundary_groups_count());
    for(const std::string& bound : mesh.boundary_names()) {
        auto& coords = quad_coords[bound];
        coords.resize(quad_shifts.at(bound).back(), std::array<T, 2>{});
        for(size_t e = 0; e < mesh.elements_count(bound); ++e) {
            const auto& el = mesh.element_1d(bound, e);
            for(size_t q = 0; q < el->qnodes_count(); ++q)
                for(size_t i = 0; i < el->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 2; ++comp)
                        coords[quad_shifts.at(bound)[e]+q][comp] += mesh.node(mesh.node_number(bound, e, i))[comp] * el->qN(i, q);
        }
    }
    return quad_coords;
}

template<class T, class I>
std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> mesh_proxy<T, I>::approx_all_jacobi_matrices_bound(const mesh_2d<T, I>& mesh,
                                                                                                                       const std::unordered_map<std::string_view, std::vector<I>>& quad_shifts) {
    check_shifts(mesh, quad_shifts);
    std::unordered_map<std::string_view, std::vector<std::array<T, 2>>> jacobi_matrices(mesh.boundary_groups_count());
    for(const std::string& bound : mesh.boundary_names()) {
        auto& matrices = jacobi_matrices[bound];
        matrices.resize(quad_shifts.at(bound).back(), std::array<T, 2>{});
        for(size_t e = 0; e < mesh.elements_count(bound); ++e) {
            const auto& el = mesh.element_1d(bound, e);
            for(size_t q = 0; q < el->qnodes_count(); ++q)
                for(size_t i = 0; i < el->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 2; ++comp)
                        matrices[quad_shifts.at(bound)[e]+q][comp] += mesh.node(mesh.node_number(bound, e, i))[comp] * el->qNxi(i, q);
        }
    }
    return jacobi_matrices;
}

template<class T, class I>
void mesh_proxy<T, I>::first_last_node_init() {
    for(size_t rank = 0, size = MPI_utils::MPI_size(); rank < size; ++rank)
        _ranges.range(rank) = { mesh().nodes_count() / size *  rank,
                                mesh().nodes_count() / size * (rank+1) + (rank == size-1) * mesh().nodes_count() % size };
}

template<class T, class I>
std::vector<std::array<T, 2>> mesh_proxy<T, I>::approx_centres_of_elements(const mesh_2d<T, I>& mesh) {
    std::vector<std::array<T, 2>> centres(mesh.elements_count(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(mesh, centres)
    for(size_t e = 0; e < centres.size(); ++e) {
        const auto& el = mesh.element_2d(e);
        const T x0 = is_trinagle(mesh.element_2d_type(e)) ? T{1}/T{3} : T{0};
        for(size_t node = 0; node < el->nodes_count(); ++node) {
            using namespace metamath::function;
            centres[e] += mesh.node(mesh.node_number(e, node)) * el->N(node, {x0, x0});
        }
    }
    return centres;
}

template<class T, class I>
void mesh_proxy<T, I>::find_elements_neighbors(std::vector<std::vector<I>>& elements_neighbors,
                                               const mesh_2d<T, I>& mesh,
                                               const std::vector<std::array<T, 2>>& centres,
                                               const std::unordered_set<I>& elements, const T r) {
    elements_neighbors.resize(mesh.elements_count());
    for(const I eL : elements)
        if (elements_neighbors[eL].empty()) {
            elements_neighbors[eL].reserve(mesh.elements_count());
            for(size_t eNL = 0; eNL < mesh.elements_count(); ++eNL)
                if(utils::distance(centres[eL], centres[eNL]) < r)
                    elements_neighbors[eL].push_back(eNL);
            elements_neighbors[eL].shrink_to_fit();
        }
}

template<class T, class I>
void mesh_proxy<T, I>::find_neighbours(const T r, const balancing_t balancing) {
    _elements_neighbors.clear();
    const std::vector<std::array<T, 2>> centres = approx_centres_of_elements(*_mesh);
    std::unordered_set<I> elements;
    for(size_t node = first_node(); node < last_node(); ++node)
        for(const I e : nodes_elements_map(node))
            elements.insert(e);
    find_elements_neighbors(_elements_neighbors, *_mesh, centres, elements, r);

    if (MPI_utils::MPI_size() > 1 && balancing != balancing_t::NO) {
        std::vector<int> data_count_per_nodes(mesh().nodes_count(), 0);
        switch(balancing) {
            case balancing_t::MEMORY: {
                std::vector<bool> nodes_flags(mesh().nodes_count());
#pragma omp parallel for default(none) shared(data_count_per_nodes) firstprivate(nodes_flags) schedule(dynamic)
                for(size_t node = first_node(); node < last_node(); ++node) {
                    std::fill(std::next(nodes_flags.begin(), node), nodes_flags.end(), false);
                    for(const I eL : nodes_elements_map(node))
                        for(const I eNL : neighbors(eL))
                            for(size_t j = 0; j < mesh().element_2d(eNL)->nodes_count(); ++j)
                                if (node <= mesh().node_number(eNL, j))
                                    if (!nodes_flags[mesh().node_number(eNL, j)]) {
                                        ++data_count_per_nodes[node];
                                        nodes_flags[mesh().node_number(eNL, j)] = true;
                                    }
                }
            }
            break;

            case balancing_t::SPEED:
#pragma omp parallel for default(none) shared(data_count_per_nodes) schedule(dynamic)
                for(size_t node = first_node(); node < last_node(); ++node)
                    for(const I eL : nodes_elements_map(node))
                        for(const I eNL : neighbors(eL))
                            for(size_t j = 0; j < mesh().element_2d(eNL)->nodes_count(); ++j)
                                if (node <= mesh().node_number(eNL, j))
                                    ++data_count_per_nodes[node];
            break;

            default:
                throw std::invalid_argument{"Unknown balancing type"};
        }

        data_count_per_nodes = MPI_utils::all_to_all<int>(data_count_per_nodes, _ranges);

        size_t sum = 0, curr_rank = 0;
        const size_t mean = std::accumulate(data_count_per_nodes.cbegin(), data_count_per_nodes.cend(), size_t{0}) / MPI_utils::MPI_size();
        for(size_t node = 0; node < data_count_per_nodes.size(); ++node) {
            sum += data_count_per_nodes[node];
            if (sum > mean) {
                _ranges.range(curr_rank).back() = node;
                _ranges.range(curr_rank+1).front() = node;
                ++curr_rank;
                sum = 0;
            }
        }

        elements.clear();
        for(size_t node = first_node(); node < last_node(); ++node)
            for(const I e : nodes_elements_map(node))
                elements.insert(e);
        find_elements_neighbors(_elements_neighbors, *_mesh, centres, elements, r);
        for(size_t e = 0; e < _elements_neighbors.size(); ++e)
            if (elements.find(e) == elements.cend()) {
                _elements_neighbors[e].clear();
                _elements_neighbors[e].shrink_to_fit();
            }
    }
}

template<class T, class I>
template<class Vector>
T mesh_proxy<T, I>::integrate_solution(const Vector& sol) const {
    if(mesh().nodes_count() != size_t(sol.size()))
        throw std::logic_error{"mesh.nodes_count() != T.size()"};

    T integral = 0;
    for(size_t e = 0; e < mesh().elements_count(); ++e) {
        const auto& el = mesh().element_2d(e);
        auto J = jacobi_matrix(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J)
            for(size_t i = 0; i < el->nodes_count(); ++i) {
                integral += el->weight(q) * el->qN(i, q) * sol[mesh().node_number(e, i)] * jacobian(*J);
        }
    }
    return integral;
}

template<class T, class I>
template<class Vector>
std::array<std::vector<T>, 2> mesh_proxy<T, I>::gradient(const Vector& sol) const {
    if(mesh().nodes_count() != sol.size())
        throw std::logic_error{"mesh.nodes_count() != sol.size()"};

    std::vector<T> dx(sol.size(), 0), dy(sol.size(), 0);
#pragma omp parallel for default(none) shared(dx, dy, sol)
    for(size_t node = 0; node < mesh().nodes_count(); ++node) {
        T node_area = T{0};
        for(const I e : nodes_elements_map(node)) {
            const auto& el = mesh().element_2d(e);
            const size_t i = global_to_local_numbering(e).find(node)->second;
            const std::array<T, 4>& J = jacobi_matrix_node(node);
            const T jac = jacobian(J);
            for(size_t j = 0; j < el->nodes_count(); ++j) {
                const T Nxi  = el->Nxi (j, el->node(i)) * element_area(e), // TODO: optimize Nxi and Neta calculations
                        Neta = el->Neta(j, el->node(i)) * element_area(e);
                dx[node] += ( J[3] * Nxi - J[2] * Neta) * sol[mesh().node_number(e, j)] / jac;
                dy[node] += (-J[1] * Nxi + J[0] * Neta) * sol[mesh().node_number(e, j)] / jac;
            }
            node_area += element_area(e);
        }
        dx[node] /= node_area;
        dy[node] /= node_area;
    }
    return {std::move(dx), std::move(dy)};
}

template<class T, class I>
std::vector<T> mesh_proxy<T, I>::approx_in_quad(const std::vector<T>& x) {
    std::vector<T> x_in_quad(quad_shift(mesh().elements_count()), 0);
#pragma omp parallel for default(none) shared(x, x_in_quad)
    for(size_t e = 0; e < mesh().elements_count(); ++e) {
        const auto& el = mesh().element_2d(e);
        for(size_t q = 0, shift = quad_shift(e); q < el->qnodes_count(); ++q, ++shift)
            for(size_t i = 0; i < el->nodes_count(); ++i)
                x_in_quad[shift] += x[mesh().node_number(e, i)] * el->qN(i, q);
    }
    return x_in_quad;
}

}

#endif