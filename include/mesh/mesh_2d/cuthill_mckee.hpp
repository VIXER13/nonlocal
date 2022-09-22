#ifndef NONLOCAL_CUTHILL_MCKEE_HPP
#define NONLOCAL_CUTHILL_MCKEE_HPP

#include "mesh_proxy.hpp"

#include "nonlocal_constants.hpp"

#include <map>
#include <set>

namespace nonlocal::mesh {

class _cuthill_mckee final {
    explicit _cuthill_mckee() noexcept = default;

    class initializer_base {
        std::vector<bool> _is_include;

    protected:
        explicit initializer_base(const size_t size)
            : _is_include(size, false) {}

        bool check_neighbour(const size_t node, const size_t neighbour) {
            const bool result = node != neighbour && !_is_include[neighbour];
            if (result)
                _is_include[neighbour] = true;
            return result;
        }

    public:
        void reset_include() {
            std::fill(_is_include.begin(), _is_include.end(), false);
        }
    };

    template<class I>
    class shifts_initializer final : public initializer_base {
        const size_t _size = 0;

    public:
        explicit shifts_initializer(const size_t size)
            : initializer_base{size}
            , _size{size} {}

        std::vector<I> init_data_vector() const {
            return std::vector<I>(_size + 1, 0);
        }

        void add(std::vector<I>& shifts, const size_t node, const size_t neighbour) {
            if (check_neighbour(node, neighbour))
                ++shifts[node + 1];
        }
    };

    template<class I>
    class indices_initializer final : public initializer_base {
        const std::vector<I>& _shifts;
        I _node_shift = 0;

    public:
        explicit indices_initializer(const std::vector<I>& shifts)
            : initializer_base{shifts.size() - 1}
            , _shifts{shifts} {}

        void reset_include() {
            _node_shift = 0;
            initializer_base::reset_include();
        }

        std::vector<I> init_data_vector() const {
            return std::vector<I>(_shifts.back(), 0);
        }

        void add(std::vector<I>& indices, const size_t node, const size_t neighbour) {
            if (check_neighbour(node, neighbour)) {
                indices[_shifts[node] + _node_shift] = neighbour;
                ++_node_shift;
            }
        }
    };

    template<class I>
    struct node_graph final {
        std::vector<I> shifts, indices;

        size_t neighbours_count(const size_t node) const {
            return shifts[node + 1] - shifts[node];
        }
    };

    template<class I>
    static void prepare_shifts(std::vector<I>& shifts) {
        for(size_t i = 1; i < shifts.size(); ++i)
            shifts[i] += shifts[i - 1];
    }

    template<theory_t Theory, class T, class I, class Initializer>
    static std::vector<I> init_vector(const mesh_proxy<T, I>& mesh, Initializer&& initializer) {
        std::vector<I> data_vector = initializer.init_data_vector();
#pragma omp parallel for default(none) shared(mesh, data_vector) firstprivate(initializer)
        for(size_t node = mesh.first_node(); node < mesh.last_node(); ++node) {
            initializer.reset_include();
            for(const I eL : mesh.nodes_elements_map(node)) {
                if constexpr (Theory == theory_t::LOCAL)
                    for(size_t jL = 0; jL < mesh.mesh().nodes_count(eL); ++jL)
                        initializer.add(data_vector, node, mesh.mesh().node_number(eL, jL));
                if constexpr (Theory == theory_t::NONLOCAL)
                    for(const I eNL : mesh.neighbors(eL))
                        for(size_t jNL = 0; jNL < mesh.mesh().nodes_count(eNL); ++jNL)
                            initializer.add(data_vector, node, mesh.mesh().node_number(eNL, jNL));
            }
        }
        return data_vector;
    }

    template<class T, class I>
    static node_graph<I> init_graph(const mesh_proxy<T, I>& mesh, const bool is_nonlocal) {
        _cuthill_mckee::node_graph<I> graph;
        graph.shifts = is_nonlocal ? init_vector<theory_t::NONLOCAL>(mesh, shifts_initializer<I>{mesh.mesh().nodes_count()})
                                   : init_vector<theory_t::   LOCAL>(mesh, shifts_initializer<I>{mesh.mesh().nodes_count()});
        prepare_shifts(graph.shifts);
        graph.indices = is_nonlocal ? init_vector<theory_t::NONLOCAL>(mesh, indices_initializer{graph.shifts})
                                    : init_vector<theory_t::   LOCAL>(mesh, indices_initializer{graph.shifts});
        return graph;
    }

    template<class I>
    static I node_with_minimum_neighbours(const node_graph<I>& graph) {
        I curr_node = 0, min_neighbours_count = std::numeric_limits<I>::max();
        for(size_t node = 0; node < graph.shifts.size() - 1; ++node)
            if (const I neighbours_count = graph.neighbours_count(node); neighbours_count < min_neighbours_count) {
                curr_node = node;
                min_neighbours_count = neighbours_count;
            }
        return curr_node;
    }

    template<class I>
    static std::vector<I> calculate_permutation(const node_graph<I>& graph, const I init_node) {
        std::vector<I> permutation(graph.shifts.size() - 1, I{-1});
        I curr_index = 0;
        permutation[init_node] = curr_index++;
        std::unordered_set<I> curr_layer{init_node}, next_layer;
        while (curr_index < permutation.size()) {
            next_layer.clear();
            for(const I node : curr_layer) {
                std::multimap<I, I> neighbours;
                for(I shift = graph.shifts[node]; shift < graph.shifts[node+1]; ++shift)
                    if (const I neighbour_node = graph.indices[shift]; permutation[neighbour_node] == -1) {
                        next_layer.emplace(neighbour_node);
                        neighbours.emplace(graph.neighbours_count(neighbour_node), neighbour_node);
                    }
                for(const auto [_, neighbour] : neighbours)
                    permutation[neighbour] = curr_index++;
            }
            std::swap(curr_layer, next_layer);
            for(const I i : permutation)
                std::cout << i << ' ';
            std::cout << std::endl;
        }
        return permutation;
    }

public:
    template<class T, class I>
    friend std::vector<I> cuthill_mckee(const mesh_proxy<T, I>& mesh, const bool is_nonlocal);
};

template<class T, class I>
std::vector<I> cuthill_mckee(const mesh_proxy<T, I>& mesh, const bool is_nonlocal) {
    const _cuthill_mckee::node_graph<I> graph = _cuthill_mckee::init_graph(mesh, is_nonlocal);
    return _cuthill_mckee::calculate_permutation(graph, _cuthill_mckee::node_with_minimum_neighbours(graph));
}

}

#endif