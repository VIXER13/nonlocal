#pragma once

#include <metamath/finite_elements/base/finite_element_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/geometry/side_1d.hpp>

#include <type_traits>

namespace metamath::finite_element {

template<class T>
class element_1d_base : public element_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    ~element_1d_base() override = default;

    virtual T node(const size_t i) const = 0;

    virtual T N  (const size_t i, const T xi) const = 0;
    virtual T Nxi(const size_t i, const T xi) const = 0;

    virtual T boundary(const side_1d bound) const = 0;
};

}