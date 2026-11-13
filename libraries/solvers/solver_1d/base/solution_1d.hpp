#pragma once

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_1d {

template<class T>
class solution_1d {
    const std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    const std::vector<model_parameters<1, T>> _models;

    std::vector<T> mix_nonlocal(const std::vector<T>& local) const;

  protected:
    explicit solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, std::vector<model_parameters<1, T>>&& models);

    template<class Parameter>
    std::vector<T> weighted_gradient(const std::vector<T>& values, const std::vector<Parameter>& parameters,
                                     const coefficient_t<T, 1> Parameter::* coefficient, const T sign) const;

  public:
    virtual ~solution_1d() noexcept = default;

    const mesh::mesh_1d<T>& mesh() const noexcept;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const model_parameters<1, T>& model(const size_t segment) const noexcept;
};

template<class T>
solution_1d<T>::solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, std::vector<model_parameters<1, T>>&& models)
    : _mesh{ mesh }, _models{ std::move(models) } {}

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

template<class T>
std::vector<T> solution_1d<T>::mix_nonlocal(const std::vector<T>& local) const {
    const auto& el = mesh().element();
    auto result = local;
    for (const size_t segment : mesh().segments()) {
        const auto& [influence, local_weight] = model(segment);
        if (theory_type(local_weight) == theory_t::LOCAL)
            continue;
        const T integration_weight = nonlocal_weight(local_weight) * mesh().jacobian(segment);
        for (const size_t eL : mesh().elements(segment))
            for (const size_t qL : el.qnodes()) {
                const T x = mesh().qnode_coord(eL, qL);
                T integral = T{ 0 };
                for (const size_t eNL : mesh().neighbours(eL))
                    for (const size_t qNL : el.qnodes())
                        integral += el.weight(qNL) * influence(x, mesh().qnode_coord(eNL, qNL)) *
                                    local[mesh().qnode_number(eNL, qNL)];
                const size_t index = mesh().qnode_number(eL, qL);
                result[index] = local_weight * local[index] + integration_weight * integral;
            }
    }
    return result;
}

template<class T>
template<class Parameter>
std::vector<T> solution_1d<T>::weighted_gradient(const std::vector<T>& values, const std::vector<Parameter>& parameters,
                                                const coefficient_t<T, 1> Parameter::* coefficient, const T sign) const {
    static constexpr auto depends_on_solution = [](const coefficient_t<T, 1>& value) noexcept {
        return std::holds_alternative<solution_dependency<T, 1>>(value);
    };
    auto local = mesh::utils::gradient_in_qnodes(mesh(), values);
    const auto values_in_qnodes = std::ranges::any_of(parameters, depends_on_solution, coefficient)
                                      ? mesh::utils::from_nodes_to_qnodes(mesh(), values)
                                      : std::vector<T>{};
    for (const size_t segment : mesh().segments()) {
        const auto& material_coefficient = parameters[segment].*coefficient;
        const bool is_solution_dependent = depends_on_solution(material_coefficient);
        for (const size_t e : mesh().elements(segment))
            for (const size_t q : mesh().element().qnodes()) {
                const size_t index = mesh().qnode_number(e, q);
                const point<T, 1> x{ mesh().qnode_coord(e, q) };
                const T value = is_solution_dependent ? values_in_qnodes[index] : T{ 0 };
                local[index] *= sign * nonlocal::evaluate(material_coefficient, x, value);
            }
    }
    return mix_nonlocal(local);
}

} // namespace nonlocal::solver_1d
