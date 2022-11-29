#ifndef NONLOCAL_ELEMENTS_SET_2D_HPP
#define NONLOCAL_ELEMENTS_SET_2D_HPP

#include "metamath.hpp"

namespace nonlocal::mesh {

template<class T>
using element_integrate_1d = metamath::finite_element::element_1d_integrate_base<T>;
template<class T>
using element_integrate_2d = metamath::finite_element::element_2d_integrate_base<T>;

template<class T>
using finite_element_1d_sptr = std::shared_ptr<element_integrate_1d<T>>;
template<class T>
using finite_element_2d_sptr = std::shared_ptr<element_integrate_2d<T>>;

template<class T>
class elements_set {
protected:
    std::vector<finite_element_1d_sptr<T>> _elements_1d;
    std::vector<finite_element_2d_sptr<T>> _elements_2d;
    const std::unordered_map<size_t, size_t> _model_to_local_1d;
    const std::unordered_map<size_t, size_t> _model_to_local_2d;
    const std::vector<size_t> _local_to_model_1d;
    const std::vector<size_t> _local_to_model_2d;

    static std::vector<size_t> local_to_model(const std::unordered_map<size_t, size_t>& model_to_local) {
        std::vector<size_t> result(model_to_local.size());
        for(const auto& [model, local] : model_to_local)
            result[local] = model;
        return result;
    }

    explicit elements_set(std::vector<finite_element_1d_sptr<T>>&& elements_1d,
                          std::vector<finite_element_2d_sptr<T>>&& elements_2d,
                          std::unordered_map<size_t, size_t>&& model_to_local_1d,
                          std::unordered_map<size_t, size_t>&& model_to_local_2d)
        : _elements_1d{std::move(elements_1d)}
        , _elements_2d{std::move(elements_2d)}
        , _model_to_local_1d{std::move(model_to_local_1d)}
        , _model_to_local_2d{std::move(model_to_local_2d)}
        , _local_to_model_1d{local_to_model(_model_to_local_1d)}
        , _local_to_model_2d{local_to_model(_model_to_local_2d)} {}

public:
    virtual ~elements_set() noexcept = default;

    element_integrate_1d<T>& element_1d(const size_t local) {
        return *_elements_1d[local];
    }

    const element_integrate_1d<T>& element_1d(const size_t local) const {
        return *_elements_1d[local];
    }

    element_integrate_2d<T>& element_2d(const size_t local) {
        return *_elements_2d[local];
    }

    const element_integrate_2d<T>& element_2d(const size_t local) const {
        return *_elements_2d[local];
    }

    size_t model_to_local_1d(const size_t model) const {
        return _model_to_local_1d.at(model);
    }

    size_t model_to_local_2d(const size_t model) const {
        return _model_to_local_2d.at(model);
    }

    size_t local_to_model_1d(const size_t local) const {
        return _local_to_model_1d.at(local);
    }

    size_t local_to_model_2d(const size_t local) const {
        return _local_to_model_2d.at(local);
    }
};

}

#endif