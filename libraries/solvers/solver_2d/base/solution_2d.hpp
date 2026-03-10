#pragma once

#include <mesh/mesh_2d/mesh_2d.hpp>
#include <constants/nonlocal_constants.hpp>
#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d {

template<class T, class I>
class solution_2d {
    const std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    const std::unordered_map<std::string, model_parameters<2, T>> _models;

protected:
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                         const std::unordered_map<std::string, model_parameters<2, T>>& models);

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

}