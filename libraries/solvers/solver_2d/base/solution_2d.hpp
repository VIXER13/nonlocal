#pragma once

#include <mesh/mesh_2d/mesh_2d.hpp>
#include <constants/nonlocal_constants.hpp>
#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d {

template<class T>
class solution_2d {
    const std::shared_ptr<mesh::mesh_2d<T>> _mesh;
    const std::unordered_map<std::string, model_parameters<2, T>> _models;

protected:
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                         const std::unordered_map<std::string, model_parameters<2, T>>& models);

public:
    virtual ~solution_2d() noexcept = default;

    const mesh::mesh_2d<T>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T>>& mesh_ptr() const noexcept;
    const model_parameters<2, T>& model(const std::string& group) const;
};

template<class T>
solution_2d<T>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _mesh{mesh} {}

template<class T>
solution_2d<T>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                               const std::unordered_map<std::string, model_parameters<2, T>>& models)
    : _mesh{mesh}
    , _models{models} {}

template<class T>
const mesh::mesh_2d<T>& solution_2d<T>::mesh() const {
    return *_mesh;
}

template<class T>
const std::shared_ptr<mesh::mesh_2d<T>>& solution_2d<T>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T>
const model_parameters<2, T>& solution_2d<T>::model(const std::string& group) const {
    return _models.at(group);
}

}