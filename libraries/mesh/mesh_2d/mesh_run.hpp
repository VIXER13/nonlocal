#ifndef NONLOCAL_MESH_RUN_HPP
#define NONLOCAL_MESH_RUN_HPP

#include "indexator_base.hpp"
#include "mesh_2d.hpp"

#include <iterator>

namespace nonlocal::mesh {

template<std::floating_point T, std::integral I, class Nodes, class Runner>
void mesh_run(const mesh_2d<T, I>& mesh,
              const Nodes& nodes, 
              const std::unordered_map<std::string, theory_t>& theories, 
              Runner&& runner) {
#pragma omp parallel for default(none) shared(mesh, nodes, theories) firstprivate(runner) schedule(dynamic)
    for(size_t i = 0; i < nodes.size(); ++i) {
        const size_t node = nodes[i];
        if constexpr (std::is_base_of_v<indexator_base, Runner>)
            runner.reset(node);
        for(const I eL : mesh.elements(node)) {
            const size_t iL = mesh.global_to_local(eL, node);
            const std::string& group = mesh.container().group(eL);
            if (const theory_t theory = theories.at(group); theory == theory_t::LOCAL)
                for(const size_t jL : std::ranges::iota_view{0u, mesh.container().nodes_count(eL)})
                    runner(group, eL, iL, jL);
            else if (theory == theory_t::NONLOCAL)
                for(const I eNL : mesh.neighbours(eL))
                    for(const size_t jNL : std::ranges::iota_view{0u, mesh.container().nodes_count(eNL)})
                        runner(group, eL, eNL, iL, jNL);
            else
                throw std::domain_error{"Unknown theory."};
        }
    }
}

}

#endif