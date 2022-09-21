#ifndef NONLOCAL_CUTHILL_MCKEE_HPP
#define NONLOCAL_CUTHILL_MCKEE_HPP

#include "mesh_proxy.hpp"

#include <map>
#include <set>

namespace nonlocal::mesh {

class _cuthill_mckee final {
    explicit _cuthill_mckee() noexcept = default;

    enum class add_t : bool {SHIFT, NEIGHBOUR};

    template<class I>
    struct node_graph final {
        std::vector<I> shifts, indices;
        std::vector<bool> is_include;
        size_t node_shift = 0;

        explicit node_graph(const size_t size)
            : shifts(size + 1, I{0})
            , is_include(size, false) {}

        void reset_include() {
            node_shift = 0;
            std::fill(is_include.begin(), is_include.end(), false);
        }

        bool check_neighbour(const size_t node, const size_t neighbour) {
            const bool result = node != neighbour && !is_include[neighbour];
            if (result)
                is_include[neighbour] = true;
            return result;
        }

        template<add_t Add>
        void add(const size_t node, const size_t neighbour) {
            if (check_neighbour(node, neighbour)) {
                if constexpr (Add == add_t::SHIFT)
                    ++shifts[node+1];
                if constexpr (Add == add_t::NEIGHBOUR) {
                    indices[shifts[node] + node_shift] = neighbour;
                    ++node_shift;
                }
            }
        }

        void prepare_memory() {
            for(size_t i = 1; i < shifts.size(); ++i)
                shifts[i] += shifts[i - 1];
            indices.resize(shifts.back(), I{-1});
        }

        size_t neighbours_count(const size_t node) const {
            return shifts[node + 1] - shifts[node];
        }
    };

    template<bool is_local, add_t Add, class T, class I>
    static void mesh_run(node_graph<I>& graph, const mesh_proxy<T, I>& mesh) {
        for(size_t node = mesh.first_node(); node < mesh.last_node(); ++node) {
            graph.reset_include();
            for(const I eL : mesh.nodes_elements_map(node)) {
                if constexpr (is_local)
                    for(size_t jL = 0; jL < mesh.mesh().nodes_count(eL); ++jL)
                        graph.template add<Add>(node, mesh.mesh().node_number(eL, jL));
                if constexpr (!is_local)
                    for(const I eNL : mesh.neighbors(eL))
                        for(size_t jNL = 0; jNL < mesh.mesh().nodes_count(eNL); ++jNL)
                            graph.template add<Add>(node, mesh.mesh().node_number(eNL, jNL));
            }
        }
    }

    template<class I>
    static I find_node_with_minimum_neighbours(const node_graph<I>& graph) {
        I curr_node = 0, min_neighbours_count = std::numeric_limits<I>::max();
        for(size_t node = 0; node < graph.shifts.size() - 1; ++node)
            if (const I neighbours_count = graph.neighbours_count(node); neighbours_count < min_neighbours_count) {
                curr_node = node;
                min_neighbours_count = neighbours_count;
            }
        return curr_node;
    }

    template<class T, class I>
    static node_graph<I> init_graph(const mesh_proxy<T, I>& mesh, const bool is_nonlocal) {
        _cuthill_mckee::node_graph<I> graph(mesh.mesh().nodes_count());
        if (is_nonlocal) _cuthill_mckee::mesh_run<false, _cuthill_mckee::add_t::SHIFT>(graph, mesh);
        else             _cuthill_mckee::mesh_run<true , _cuthill_mckee::add_t::SHIFT>(graph, mesh);
        graph.prepare_memory();
        if (is_nonlocal) _cuthill_mckee::mesh_run<false, _cuthill_mckee::add_t::NEIGHBOUR>(graph, mesh);
        else             _cuthill_mckee::mesh_run<true , _cuthill_mckee::add_t::NEIGHBOUR>(graph, mesh);
        return graph;
    }

public:
    template<class T, class I>
    friend std::vector<I> cuthill_mckee(const mesh_proxy<T, I>& mesh, const bool is_nonlocal);
};

template<class T, class I>
std::vector<I> cuthill_mckee(const mesh_proxy<T, I>& mesh, const bool is_nonlocal) {
    const _cuthill_mckee::node_graph<I> graph = _cuthill_mckee::init_graph(mesh, is_nonlocal);
    const I init_node = _cuthill_mckee::find_node_with_minimum_neighbours(graph);

    std::vector<I> permutation(mesh.mesh().nodes_count(), I{-1});
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
    }
    return permutation;
}

}

#endif