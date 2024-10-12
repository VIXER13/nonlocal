#pragma once

#include "mesh_2d.hpp"
#include "nonlocal_constants.hpp"
#include "../equation_parameters.hpp"

namespace nonlocal {

template<class T, class I>
class solution_2d {
    const std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    const std::unordered_map<std::string, model_parameters<2, T>> _models;

protected:
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                         const std::unordered_map<std::string, model_parameters<2, T>>& models);

    // template<class Callback>
    // void calc_nonlocal(const Callback& callback) const;

public:
    virtual ~solution_2d() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    const model_parameters<2, T>& model(const std::string& group) const;
};

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                               const std::unordered_map<std::string, model_parameters<2, T>>& models)
    : _mesh{mesh}
    , _models{models} {}

template<class T, class I>
const mesh::mesh_2d<T, I>& solution_2d<T, I>::mesh() const {
    return *_mesh;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_2d<T, I>>& solution_2d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
const model_parameters<2, T>& solution_2d<T, I>::model(const std::string& group) const {
    return _models.at(group);
}

// template<class T, class I>
// template<class Callback>
// void solution_2d<T, I>::calc_nonlocal(const Callback& callback) const {
//     const auto process_nodes = mesh().process_nodes();
//     std::vector<bool> neighbors(mesh().container().elements_count(), true);
// #pragma omp parallel for default(none) shared(process_nodes) firstprivate(neighbors, callback)
//     for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
//         std::fill(neighbors.begin(), neighbors.end(), true);
//         for(const I eL : mesh().elements(node))
//             for(const I eNL : mesh().neighbours(eL))
//                 if (neighbors[eNL]) {
//                     callback(eNL, node);
//                     neighbors[eNL] = false;
//                 }
//     }
// }

}