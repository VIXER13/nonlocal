#pragma once

#include "mesh_2d_utils.hpp"

namespace nonlocal::mesh::utils {

class _cuthill_mckee final {
    explicit _cuthill_mckee() noexcept = default;

    class neighbour_checker : public indexator_base {
        std::vector<bool> _included;

    protected:
        explicit neighbour_checker(const size_t size)
            : indexator_base{false}
            , _included(size, false) {}

        bool check_neighbour(const size_t node, const size_t neighbour) {
            const bool result = node != neighbour && !_included[neighbour];
            if (result)
                _included[neighbour] = true;
            return result;
        }

    public:
        void reset(const size_t) override {
            std::fill(_included.begin(), _included.end(), false);
        }
    };

    template<class T, class I>
    class shifts_initializer final : public neighbour_checker {
        std::vector<size_t>& _shifts;
        const mesh_container_2d<T, I>& _mesh;

        void check(const size_t node, const size_t neighbour) {
            if (check_neighbour(node, neighbour))
                ++_shifts[node + 1];
        }

    public:
        explicit shifts_initializer(std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh)
            : neighbour_checker{shifts.size() - 1}
            , _shifts{shifts}
            , _mesh{mesh} {}

        void operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
            check(_mesh.node_number(e, i), _mesh.node_number(e, j));
        }

        void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            check(_mesh.node_number(eL, iL), _mesh.node_number(eNL, jNL));
        }
    };

    template<class T, class I>
    class indices_initializer final : public neighbour_checker {
        const std::vector<size_t>& _shifts;
        std::vector<I>& _indices;
        I _node_shift = 0;
        const mesh_container_2d<T, I>& _mesh;

        void check(const size_t node, const size_t neighbour) {
            if (check_neighbour(node, neighbour)) {
                _indices[_shifts[node] + _node_shift] = neighbour;
                ++_node_shift;
            }
        }

    public:
        explicit indices_initializer(std::vector<I>& indices, const std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh)
            : neighbour_checker{shifts.size() - 1}
            , _shifts{shifts}
            , _indices{indices}
            , _mesh{mesh} {}

        void reset(const size_t) override {
            _node_shift = 0;
            neighbour_checker::reset(0);
        }

        void operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
            check(_mesh.node_number(e, i), _mesh.node_number(e, j));
        }

        void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            check(_mesh.node_number(eL, iL), _mesh.node_number(eNL, jNL));
        }
    };

    template<class I>
    struct node_graph final {
        std::vector<size_t> shifts;
        std::vector<I> indices;

        size_t neighbours_count(const size_t node) const {
            return shifts[node + 1] - shifts[node];
        }
    };

    template<class I>
    static void accumulate_shifts(std::vector<I>& shifts) {
        for(const size_t i : std::ranges::iota_view{0u, shifts.size()})
            shifts[i] += shifts[i - 1];
    }

    template<class T, class I>
    static node_graph<I> init_graph(const mesh_2d<T, I>& mesh, const std::unordered_map<std::string, theory_t>& theories) {
        _cuthill_mckee::node_graph<I> graph;
        graph.shifts.resize(mesh.container().nodes_count() + 1, 0);
        mesh_run(mesh, mesh.container().nodes(), theories, shifts_initializer{graph.shifts, mesh.container()});
        accumulate_shifts(graph.shifts);
        graph.indices.resize(graph.shifts.back());
        mesh_run(mesh, mesh.container().nodes(), theories, indices_initializer{graph.indices, graph.shifts, mesh.container()});
        return graph;
    }

    template<class I>
    static I node_with_minimum_neighbours(const node_graph<I>& graph) {
        I curr_node = 0, min_neighbours_count = std::numeric_limits<I>::max();
        for(const size_t node : std::ranges::iota_view{0u, graph.shifts.size() - 1})
            if (const I neighbours_count = graph.neighbours_count(node); neighbours_count < min_neighbours_count) {
                curr_node = node;
                min_neighbours_count = neighbours_count;
            }
        return curr_node;
    }

    template<class I>
    static std::vector<size_t> calculate_permutation(const node_graph<I>& graph, const I init_node) {
        std::vector<size_t> permutation(graph.shifts.size() - 1, I{-1});
        I curr_index = 0;
        permutation[init_node] = curr_index++;
        std::unordered_set<I> curr_layer{init_node}, next_layer;
        while (curr_index < permutation.size()) {
            next_layer.clear();
            for(const I node : curr_layer) {
                std::multimap<I, I> neighbours;
                for(const I shift : std::ranges::iota_view{graph.shifts[node], graph.shifts[node + 1]})
                    if (const I neighbour_node = graph.indices[shift]; permutation[neighbour_node] == I{-1})
                        neighbours.emplace(graph.neighbours_count(neighbour_node), neighbour_node);
                for(const auto [_, neighbour] : neighbours) {
                    next_layer.emplace(neighbour);
                    permutation[neighbour] = curr_index++;
                }
            }
            std::swap(curr_layer, next_layer);
        }
        return permutation;
    }

public:
    template<class T, class I>
    friend std::vector<size_t> cuthill_mckee(const mesh_2d<T, I>& mesh, const bool only_local, const bool reverse);
};

template<class I>
std::vector<I> reverse_permutation(const std::vector<I>& permutation) {
    std::vector<I> reversed(permutation.size());
    I index = 0;
    for(const I i : permutation)
        reversed[i] = index++;
    return reversed;
}

template<class T, class I>
std::vector<size_t> cuthill_mckee(const mesh_2d<T, I>& mesh, const bool only_local, const bool reverse) {
    const _cuthill_mckee::node_graph<I> graph = _cuthill_mckee::init_graph(mesh, theories(mesh, only_local));
    std::vector<size_t> permutation = _cuthill_mckee::calculate_permutation(graph, _cuthill_mckee::node_with_minimum_neighbours(graph));
    return reverse ? reverse_permutation(permutation) : permutation;
}

}