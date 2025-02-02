#pragma once

#include <metamath/metamath.hpp>

namespace nonlocal::mesh {

template<class T>
using element_integrate_1d = metamath::finite_element::element_1d_integrate_base<T>;
template<class T>
using element_integrate_2d = metamath::finite_element::element_2d_integrate_base<T>;

template<class T>
using finite_element_1d_sptr = std::shared_ptr<element_integrate_1d<T>>;
template<class T>
using finite_element_2d_sptr = std::shared_ptr<element_integrate_2d<T>>;

enum class element_1d_t : uint8_t {
    LINEAR,
    QUADRATIC
};

enum class element_2d_t : uint8_t {
    TRIANGLE,
    QUADRATIC_TRIANGLE,
    BILINEAR,
    QUADRATIC_SERENDIPITY,
    QUADRATIC_LAGRANGE
};

template<class T>
class elements_set {
protected:
    std::vector<finite_element_1d_sptr<T>> _elements_1d;
    std::vector<finite_element_2d_sptr<T>> _elements_2d;
    const std::unordered_map<size_t, element_1d_t> _model_to_local_1d;
    const std::unordered_map<size_t, element_2d_t> _model_to_local_2d;
    const std::vector<size_t> _local_to_model_1d;
    const std::vector<size_t> _local_to_model_2d;

    template<class I>
    static std::vector<size_t> local_to_model(const std::unordered_map<size_t, I>& model_to_local) {
        std::vector<size_t> result(model_to_local.size());
        for(const auto& [model, local] : model_to_local)
            result[size_t(local)] = model;
        return result;
    }

    explicit elements_set(std::vector<finite_element_1d_sptr<T>>&& elements_1d,
                          std::vector<finite_element_2d_sptr<T>>&& elements_2d,
                          std::unordered_map<size_t, element_1d_t>&& model_to_local_1d,
                          std::unordered_map<size_t, element_2d_t>&& model_to_local_2d)
        : _elements_1d{std::move(elements_1d)}
        , _elements_2d{std::move(elements_2d)}
        , _model_to_local_1d{std::move(model_to_local_1d)}
        , _model_to_local_2d{std::move(model_to_local_2d)}
        , _local_to_model_1d{local_to_model(_model_to_local_1d)}
        , _local_to_model_2d{local_to_model(_model_to_local_2d)} {}

public:
    virtual ~elements_set() noexcept = default;

    element_integrate_1d<T>& element_1d(const element_1d_t local) {
        return *_elements_1d[size_t(local)];
    }

    const element_integrate_1d<T>& element_1d(const element_1d_t local) const {
        return *_elements_1d[size_t(local)];
    }

    element_integrate_2d<T>& element_2d(const element_2d_t local) {
        return *_elements_2d[size_t(local)];
    }

    const element_integrate_2d<T>& element_2d(const element_2d_t local) const {
        return *_elements_2d[size_t(local)];
    }

    element_1d_t model_to_local_1d(const size_t model) const {
        return _model_to_local_1d.at(model);
    }

    element_2d_t model_to_local_2d(const size_t model) const {
        return _model_to_local_2d.at(model);
    }

    size_t local_to_model_1d(const element_1d_t local) const {
        return _local_to_model_1d[size_t(local)];
    }

    size_t local_to_model_2d(const element_2d_t local) const {
        return _local_to_model_2d[size_t(local)];
    }

    bool is_element_1d(const size_t model) const {
        return _model_to_local_1d.contains(model);
    }

    bool is_element_2d(const size_t model) const {
        return _model_to_local_2d.contains(model);
    }
};

}