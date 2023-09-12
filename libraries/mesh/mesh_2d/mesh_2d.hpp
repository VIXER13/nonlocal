#ifndef NONLOCAL_MESH_2D_HPP
#define NONLOCAL_MESH_2D_HPP

#include "mesh_container_2d_utils.hpp"

#include "MPI_utils.hpp"

namespace nonlocal::mesh {

enum class balancing_t : uint8_t { NO, MEMORY, SPEED };

template<class T>
constexpr T jacobian(const std::array<T, 2>& J) {
    return std::sqrt(J[X] * J[X] + J[Y] * J[Y]);
}

template<class T>
constexpr T jacobian(const metamath::types::square_matrix<T, 2>& J) noexcept {
    return std::abs(J[X][X] * J[Y][Y] - J[X][Y] * J[Y][X]);
}

template<class T, class I>
class mesh_2d final {
    mesh_container_2d<T, I> _mesh;

    std::vector<std::vector<I>> _node_elements;
    std::vector<std::unordered_map<I, uint8_t>> _global_to_local;

    std::vector<I> _quad_shifts;
    std::vector<std::array<T, 2>> _quad_coords;
    std::vector<metamath::types::square_matrix<T, 2>> _jacobi_matrices;

    std::vector<I> _quad_node_shift;
    std::vector<std::array<T, 2>> _derivatives;

    parallel_utils::MPI_ranges _MPI_ranges;

    std::vector<std::vector<I>> _elements_neighbors;

    T area(const std::ranges::iota_view<size_t, size_t> elements) const;

public:
    explicit mesh_2d(const std::filesystem::path& path_to_mesh);

    const mesh_container_2d<T, I>& container() const;
    
    const std::vector<I>& elements(const size_t node) const;
    size_t global_to_local(const size_t e, const size_t node) const;

    size_t quad_shift(const size_t e) const;
    std::ranges::iota_view<size_t, size_t> quad_shifts_count(const size_t e) const;
    const std::array<T, 2>& quad_coord(const size_t qshift) const;
    const std::array<T, 2>& quad_coord(const size_t e, const size_t q) const;
    const metamath::types::square_matrix<T, 2>& jacobi_matrix(const size_t qshift) const;
    const metamath::types::square_matrix<T, 2>& jacobi_matrix(const size_t e, const size_t q) const;

    size_t quad_node_shift(const size_t e, const size_t i) const;
    const std::array<T, 2>& derivatives(const size_t qshift) const;
    const std::array<T, 2>& derivatives(const size_t qnode_shift, const size_t q) const;
    const std::array<T, 2>& derivatives(const size_t e, const size_t i, const size_t q) const;

    const parallel_utils::MPI_ranges& MPI_ranges() const noexcept;
    std::ranges::iota_view<size_t, size_t> process_nodes(const size_t process = parallel_utils::MPI_rank()) const;
    std::unordered_set<I> process_elements(const size_t process = parallel_utils::MPI_rank()) const;

    const std::vector<I>& neighbours(const size_t e) const;

    T area(const size_t e) const;
    T area(const std::string& element_group) const;
    T area() const;

    void find_neighbours(const std::unordered_map<std::string, T>& radii, const balancing_t balancing = balancing_t::MEMORY, const bool add_diam = true);

    void clear();
};

template<class T, class I>
mesh_2d<T, I>::mesh_2d(const std::filesystem::path& path_to_mesh)
    : _mesh{path_to_mesh}
    , _node_elements{utils::node_elements_2d(container())}
    , _global_to_local{utils::global_to_local(container())}
    , _quad_shifts{utils::elements_quadrature_shifts_2d(container())}
    , _quad_coords{utils::approx_all_quad_nodes(container(), _quad_shifts)}
    , _jacobi_matrices{utils::approx_all_jacobi_matrices(container(), _quad_shifts)}
    , _quad_node_shift{utils::element_node_shits_quadrature_shifts_2d(container())}
    , _derivatives{utils::derivatives_in_quad(container(), _quad_shifts, _quad_node_shift, _jacobi_matrices)}
    , _MPI_ranges{container().nodes_count()}
    , _elements_neighbors(container().elements_2d_count()) {}

template<class T, class I>
const mesh_container_2d<T, I>& mesh_2d<T, I>::container() const {
    return _mesh;
}

template<class T, class I>
const std::vector<I>& mesh_2d<T, I>::elements(const size_t node) const {
    return _node_elements[node];
}

template<class T, class I>
size_t mesh_2d<T, I>::global_to_local(const size_t e, const size_t node) const {
    return _global_to_local[e].at(node);
}

template<class T, class I>
size_t mesh_2d<T, I>::quad_shift(const size_t e) const {
    return _quad_shifts[e];
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_2d<T, I>::quad_shifts_count(const size_t e) const {
    return {quad_shift(e), quad_shift(e + 1)};
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::quad_coord(const size_t qshift) const {
    return _quad_coords[qshift];
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::quad_coord(const size_t e, const size_t q) const {
    return quad_coord(quad_shift(e) + q);
}

template<class T, class I>
const metamath::types::square_matrix<T, 2>& mesh_2d<T, I>::jacobi_matrix(const size_t qshift) const {
    return _jacobi_matrices[qshift];
}

template<class T, class I>
const metamath::types::square_matrix<T, 2>& mesh_2d<T, I>::jacobi_matrix(const size_t e, const size_t q) const {
    return jacobi_matrix(quad_shift(e) + q);
}

template<class T, class I>
size_t mesh_2d<T, I>::quad_node_shift(const size_t e, const size_t i) const {
    return _quad_node_shift[e] + i * container().element_2d(e).qnodes_count();
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::derivatives(const size_t qshift) const {
    return _derivatives[qshift];
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::derivatives(const size_t qnode_shift, const size_t q) const {
    return derivatives(qnode_shift + q);
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::derivatives(const size_t e, const size_t i, const size_t q) const {
    return derivatives(quad_node_shift(e, i), q);
}

template<class T, class I>
const parallel_utils::MPI_ranges& mesh_2d<T, I>::MPI_ranges() const noexcept {
    return _MPI_ranges;
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_2d<T, I>::process_nodes(const size_t process) const {
    return _MPI_ranges.get(process);
}

template<class T, class I>
std::unordered_set<I> mesh_2d<T, I>::process_elements(const size_t process) const {
    std::unordered_set<I> proc_elements;
    for(const size_t node : process_nodes(process))
        for(const I e : elements(node))
            proc_elements.insert(e);
    return proc_elements;
}

template<class T, class I>
const std::vector<I>& mesh_2d<T, I>::neighbours(const size_t e) const {
    return _elements_neighbors[e];
}

template<class T, class I>
T mesh_2d<T, I>::area(const size_t e) const {
    T area = T{0};
    const auto& el = container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
        const T factor = el.weight(q) * jacobian(jacobi_matrix(e, q));
        for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
            area += factor * el.qN(i, q);
    }
    return area;
}

template<class T, class I>
T mesh_2d<T, I>::area(const std::ranges::iota_view<size_t, size_t> elements) const {
    const auto summator = [this](const T sum, const size_t e) { return sum + area(e); };
    return std::reduce(elements.begin(), elements.end(), T{0}, summator);
}

template<class T, class I>
T mesh_2d<T, I>::area(const std::string& element_group) const {
    if (!container().groups_names_2d().contains(element_group))
        throw std::domain_error{
            "It is not possible to calculate the area of the elements group, "
            "since the elements group " + element_group + " is missing"
        };
    return area(container().elements(element_group));
}

template<class T, class I>
T mesh_2d<T, I>::area() const {
    return area(container().elements_2d());
}

template<class T, class I>
void mesh_2d<T, I>::find_neighbours(const std::unordered_map<std::string, T>& radii, const balancing_t balancing, const bool add_diam) {
    if (radii.empty())
        return;
    _elements_neighbors.resize(container().elements_2d_count());
    const std::vector<std::array<T, 2>> centers = utils::approx_centers_of_elements(container());
    for(const auto& [group, radius] : radii) {
        if (radius == T{0})
            continue;
        const auto elements_range = container().elements(group);
        for(const size_t eL : elements_range) {
            _elements_neighbors[eL].reserve(elements_range.size());
            for(const size_t eNL : elements_range)
                if (metamath::functions::distance(centers[eL], centers[eNL]) <= radius)
                    _elements_neighbors[eL].push_back(eNL);
            _elements_neighbors[eL].shrink_to_fit();
        }
    }
}

template<class T, class I>
void mesh_2d<T, I>::clear() {
    _mesh.clear();
    _node_elements.clear();
    _node_elements.shrink_to_fit();
    _global_to_local.claer();
    _global_to_local.shrink_to_fit();
    _quad_shifts.clear();
    _quad_shifts.shrink_to_fit();
    _quad_coords.clear();
    _quad_coords.shrink_to_fit();
    _jacobi_matrices.clear();
    _jacobi_matrices.shrink_to_fit();
    _quad_node_shift.clear();
    _quad_node_shift.shrink_to_fit();
    _derivatives.clear();
    _derivatives.shrink_to_fit();
    _MPI_ranges = parallel_utils::MPI_ranges{0};
    _elements_neighbors.clear();
    _elements_neighbors.shrink_to_fit();
}

}

/*

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
                if(metamath::functions::distance(centres[eL], centres[eNL]) < r)
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

    if (parallel_utils::MPI_size() > 1 && balancing != balancing_t::NO) {
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

        data_count_per_nodes = parallel_utils::all_to_all<int>(data_count_per_nodes, _ranges);

        size_t sum = 0, curr_rank = 0;
        const size_t mean = std::accumulate(data_count_per_nodes.cbegin(), data_count_per_nodes.cend(), size_t{0}) / parallel_utils::MPI_size();
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

*/

#endif