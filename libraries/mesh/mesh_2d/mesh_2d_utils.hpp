#ifndef NONLOCAL_MESH_UTILS_HPP
#define NONLOCAL_MESH_UTILS_HPP

#include "nonzero_counter.hpp"
#include "integral_counter.hpp"
#include "find_neighbours.hpp"

#include <iterator>

namespace nonlocal::mesh::utils {

enum class balancing_t : uint8_t { NO, MEMORY, SPEED };

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

template<class T, class I, class Vector>
std::array<std::vector<T>, 2> gradient_in_qnodes(const mesh_2d<T, I>& mesh, const Vector& x) {
    if (mesh.container().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    const size_t quadratures_count = mesh.quad_shift(mesh.container().elements_2d_count());
    std::array<std::vector<T>, 2> gradient{
        std::vector<T>(quadratures_count, T{0}),
        std::vector<T>(quadratures_count, T{0})
    };
#pragma parallel for default(none) shared(gradient, mesh, x)
    for(size_t e = 0; e < mesh.container().elements_2d_count(); ++e) {
        const auto& el = mesh.container().element_2d(e);
        for(size_t q = 0, qshift = mesh.quad_shift(e); q < el.qnodes_count(); ++q, ++qshift) {
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()}) {
                const std::array<T, 2>& derivatives = mesh.derivatives(e, i, q);
                const T& val = x[mesh.container().node_number(e, i)];
                gradient[X][qshift] += val * derivatives[X];
                gradient[Y][qshift] += val * derivatives[Y];
            }
            const T jac = jacobian(mesh.jacobi_matrix(qshift));
            gradient[X][qshift] /= jac;
            gradient[Y][qshift] /= jac;
        }
    }
    return gradient;
}

template<class T, class I, class Vector>
std::vector<T> nodes_to_qnodes(const mesh_2d<T, I>& mesh, const Vector& x) {
    if (mesh.container().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    const size_t quadratures_count = mesh.quad_shift(mesh.container().elements_2d_count());
    std::vector<T> values(quadratures_count, T{0});
#pragma parallel for default(none) shared(gradient, mesh, x, values)
    for(size_t e = 0; e < mesh.container().elements_2d_count(); ++e) {
        const auto& el = mesh.container().element_2d(e);
        for(size_t q = 0, qshift = mesh.quad_shift(e); q < el.qnodes_count(); ++q, ++qshift) {
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()}) 
                values[qshift] += x[mesh.container().node_number(e, i)] * el.qN(i, q);
        }
    }
    return values;
}

template<class T, class I, class Vector>
std::vector<T> qnodes_to_nodes(const mesh_2d<T, I>& mesh, const Vector& x) {
    if (mesh.quad_shift(mesh.container().elements_2d_count()) != size_t(x.size()))
        throw std::logic_error{"Cannot approximate node values because vector size does not match number of quadrature nodes"};
    std::vector<T> approximation(mesh.container().nodes_count(), T{0});
#pragma parallel for default(none) shared(approximation, mesh, x)
    for(size_t node = 0; node < mesh.container().nodes_count(); ++node) {
        T node_area = T{0};
        for(const I e : mesh.elements(node)) {
            const T area = mesh.area(e);
            const size_t i = mesh.global_to_local(e, node);
            const auto& el = mesh.container().element_2d(e);
            const size_t qshift = mesh.quad_shift(e) + el.nearest_qnode(i);
            approximation[node] += area * x[qshift];
            node_area += area;
        }
        approximation[node] /= node_area;
    }
    return approximation;
}

template<class T, class I>
std::unordered_map<std::string, theory_t> theories(const mesh_2d<T, I>& mesh, const bool only_local) {
    std::unordered_map<std::string, theory_t> theories;
    for(const std::string& group : mesh.container().groups_2d())
        theories[group] = only_local             ? theory_t::LOCAL :
                          mesh.radius(group) > 0 ? theory_t::NONLOCAL :
                                                   theory_t::LOCAL;
    return theories;
}

template<class T, class I>
void balancing(mesh_2d<T, I>& mesh, const balancing_t balance, const bool only_local, const bool is_symmetric) {
    if (balance == balancing_t::NO || parallel::MPI_size() == 1)
        return;
    std::vector<size_t> nonzero_elements_count(mesh.container().nodes_count(), 0);
    if (balance == balancing_t::MEMORY)
        mesh_run(mesh, mesh.process_nodes(), theories(mesh, only_local),
            nonzero_counter{nonzero_elements_count, mesh.container(), is_symmetric});
    else if (balance == balancing_t::SPEED)
        mesh_run(mesh, mesh.process_nodes(), theories(mesh, only_local),
            integral_counter{nonzero_elements_count, mesh.container(), is_symmetric});
    else throw std::domain_error{"Unsupported balancing type"};
    nonzero_elements_count = parallel::all_to_all(nonzero_elements_count, mesh.MPI_ranges());
    mesh.MPI_ranges(parallel::uniform_ranges(nonzero_elements_count, parallel::MPI_size()));
    mesh.neighbours(find_neighbours(mesh, mesh.radii(), diam_adding::NO));
}

// template<class T, class I, std::ranges::random_access_range Vector>
// T integrate(const mesh_proxy<T, I>& mesh, const Vector& x) {
//     if (mesh.mesh().nodes_count() != size_t(x.size()))
//         throw std::logic_error{"The integral cannot be found because the vector size does not match the number of nodes."};
//     T integral = 0;
//     for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
//         auto J = mesh.jacobi_matrix(e);
//         const auto& el = mesh.mesh().element_2d(e);
//         for(size_t q = 0; q < el->qnodes_count(); ++q, ++J)
//             for(const size_t i : std::views::iota(size_t{0}, el->nodes_count()))
//                 integral += el->weight(q) * x[mesh.mesh().node_number(e, i)] * el->qN(i, q) * jacobian(*J);
//     }
//     return integral;
// }

// template<class T, class I, std::ranges::random_access_range Vector>
// std::vector<T> approximate_in_qnodes(const mesh_proxy<T, I>& mesh, const Vector& x) {
//     if (mesh.mesh().nodes_count() != size_t(x.size()))
//         throw std::logic_error{"The approximation in qnodes cannot be found because the vector size does not match the number of nodes."};
//     std::vector<T> approximation(mesh.quad_shift(mesh.mesh().elements_count()), T{0});
//     for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
//         const auto& el = mesh.mesh().element_2d(e);
//         for(size_t q = 0, qshift = mesh.quad_shift(e); q < el->qnodes_count(); ++q, ++qshift)
//             for(size_t i = 0; i < el->nodes_count(); ++i)
//                 approximation[qshift] += x[mesh.mesh().node_number(e, i)] * el->qN(i, q);
//     }
//     return approximation;
// }

}

#endif