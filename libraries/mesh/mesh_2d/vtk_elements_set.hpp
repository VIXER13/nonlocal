#pragma once

#include "elements_set.hpp"

namespace nonlocal::mesh {

enum class vtk_element_number : size_t {
    LINEAR = 3,
    QUADRATIC = 21,
    TRIANGLE = 5,
    QUADRATIC_TRIANGLE = 22,
    BILINEAR = 9,
    QUADRATIC_SERENDIPITY = 23,
    QUADRATIC_LAGRANGE = 28
};

template<std::floating_point T>
inline constexpr std::string_view vtk_data_type = std::is_same_v<T, float> ? "float" : "double";

template<class T>
class vtk_elements_set final : public elements_set<T> {
    template<class U, template<class, auto...> class Quadrature_Type, auto... Args>
    using quadrature = metamath::finite_element::quadrature_1d<U, Quadrature_Type, Args...>;
    template<class U, size_t N>
    using gauss = metamath::finite_element::gauss<U, N>;

    template<class U, size_t N>
    using lagrangian_element_1d = metamath::finite_element::lagrangian_element_1d<U, N>;
    template<class U, size_t Order>
    using element_1d = metamath::finite_element::element_1d<U, lagrangian_element_1d, Order>;

    template<class U, size_t Order>
    using triangle = metamath::finite_element::triangle<U, Order>;
    template<class U, size_t N, size_t M>
    using lagrangian_element_2d = metamath::finite_element::lagrangian_element_2d<U, N, M>;
    template<class U, template<class, auto...> class Element_Type, auto... Args>
    using element_2d = metamath::finite_element::element_2d<U, Element_Type, Args...>;

    static std::vector<element_integrate_1d<T>> make_default_1d_elements() {
        return {
            element_integrate_1d<T>{std::make_unique<element_1d<T, 1>>(), quadrature<T, gauss, 1>{}},
            element_integrate_1d<T>{std::make_unique<element_1d<T, 2>>(), quadrature<T, gauss, 2>{} }
        };
    }

    static std::vector<element_integrate_2d<T>> make_default_2d_elements() {
        return {
            element_integrate_2d<T>{ std::make_unique<element_2d<T, triangle, 1>>(), quadrature<T, gauss, 1>{}},
            element_integrate_2d<T>{std::make_unique<element_2d<T, triangle, 2>>(), quadrature<T, gauss, 2>{}},
            element_integrate_2d<T>{std::make_unique<element_2d<T, metamath::finite_element::serendipity, 1>>(), quadrature<T, gauss, 2>{}},
            element_integrate_2d<T>{std::make_unique<element_2d<T, metamath::finite_element::serendipity, 2>>(), quadrature<T, gauss, 3>{}},
            element_integrate_2d<T>{std::make_unique<element_2d<T, lagrangian_element_2d, 2, 2>>(), quadrature<T, gauss, 3>{}}
        };
    }

    static std::unordered_map<size_t, element_1d_t> vtk_to_local_1d() {
        return {
            {size_t(vtk_element_number::LINEAR),    element_1d_t::LINEAR},
            {size_t(vtk_element_number::QUADRATIC), element_1d_t::QUADRATIC}
        };
    }

    static std::unordered_map<size_t, element_2d_t> vtk_to_local_2d() {
        return {
            {size_t(vtk_element_number::TRIANGLE),              element_2d_t::TRIANGLE},
            {size_t(vtk_element_number::QUADRATIC_TRIANGLE),    element_2d_t::QUADRATIC_TRIANGLE},
            {size_t(vtk_element_number::BILINEAR),              element_2d_t::BILINEAR},
            {size_t(vtk_element_number::QUADRATIC_SERENDIPITY), element_2d_t::QUADRATIC_SERENDIPITY},
            {size_t(vtk_element_number::QUADRATIC_LAGRANGE),    element_2d_t::QUADRATIC_LAGRANGE}
        };
    }

public:
    explicit vtk_elements_set()
        : elements_set<T>{make_default_1d_elements(), 
                          make_default_2d_elements(), 
                          vtk_to_local_1d(), 
                          vtk_to_local_2d()} {}
    ~vtk_elements_set() noexcept override = default;
};

}