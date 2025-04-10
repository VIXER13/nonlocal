#pragma once

#include <mesh/mesh_1d/mesh_1d.hpp>
#include <solvers/equation_parameters.hpp>

namespace nonlocal {

template<class T>
class solution_1d {
    const std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    const std::vector<model_parameters<1, T>> _models;

protected:
    explicit solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                         std::vector<model_parameters<1, T>>&& models);

public:
    virtual ~solution_1d() noexcept = default;

    const mesh::mesh_1d<T>& mesh() const noexcept;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const model_parameters<1, T>& model(const size_t segment) const noexcept;
};

template<class T>
solution_1d<T>::solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                            std::vector<model_parameters<1, T>>&& models)
    : _mesh{mesh}
    , _models{std::move(models)} {}

template<class T>
const mesh::mesh_1d<T>& solution_1d<T>::mesh() const noexcept {
    return *_mesh;
}

template<class T>
const std::shared_ptr<mesh::mesh_1d<T>>& solution_1d<T>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T>
const model_parameters<1, T>& solution_1d<T>::model(const size_t segment) const noexcept {
    return _models[segment];
}

}