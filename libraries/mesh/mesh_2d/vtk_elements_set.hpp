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
    static std::vector<element_integrate_1d<T>> make_default_1d_elements() {
        using namespace metamath::finite_element;
        return {
            make_element_1d_integrated<T>(1, 1),
            make_element_1d_integrated<T>(2, 2)
        };
    }

    static std::vector<element_integrate_2d<T>> make_default_2d_elements() {
        using namespace metamath::finite_element;
        return {
            make_element_2d_integrated<T>(element_t::Triangle, 1, 1),
            make_element_2d_integrated<T>(element_t::Triangle, 2, 2),
            make_element_2d_integrated<T>(element_t::Serendipity, 1, 2),
            make_element_2d_integrated<T>(element_t::Serendipity, 2, 3),
            make_element_2d_integrated<T>(element_t::Lagrangian, 2, 3)
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