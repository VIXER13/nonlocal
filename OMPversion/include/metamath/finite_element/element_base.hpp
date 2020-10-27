#ifndef FINITE_ELEMENT_BASE_HPP
#define FINITE_ELEMENT_BASE_HPP

#include "element_1d_strategy/geometry_1d.hpp"
#include "element_2d_strategy/geometry_2d.hpp"

namespace metamath::finite_element {

class element_base {
public:
    virtual size_t nodes_count() const = 0; // Любой конечный элемент, вне зависимости от его размерности, имеет некоторое количество узлов.
    virtual ~element_base() noexcept = default;
};

template<class T>
class element_1d_base : public element_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    ~element_1d_base() override = default;

    virtual T node(const size_t i) const = 0;

    virtual T N  (const size_t i, const T xi) const = 0; // Обращение к i-ой функции формы в точке xi.
    virtual T Nxi(const size_t i, const T xi) const = 0; // Аналогично для производной.

    virtual T boundary(const side_1d bound) const = 0; // Геометрия элемента.
};

template<class T>
class element_2d_base : public element_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    ~element_2d_base() override = default;

    virtual const std::array<T, 2>& node(const size_t i) const = 0;

    virtual T N   (const size_t i, const std::array<T, 2>& xi) const = 0; // Обращение к i-ой функции формы в точке (xi, eta).
    virtual T Nxi (const size_t i, const std::array<T, 2>& xi) const = 0; // Аналогично для производных.
    virtual T Neta(const size_t i, const std::array<T, 2>& xi) const = 0;

    virtual T boundary(const side_2d bound, const T x) const = 0; // Геометрия элемента.
};

}

#endif