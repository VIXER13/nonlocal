#ifndef FINITE_ELEMENT_GEOMETRY_2D_HPP
#define FINITE_ELEMENT_GEOMETRY_2D_HPP

#include <array>
#include <functional>
#include "../../symdiff/symdiff.hpp"

namespace metamath::finite_element {

// Двумерную геометрию можно описать четырьями функциями, каждая из которых описывает границу интегрирования по каждой из сторон.
enum class side_2d : uint8_t { LEFT, RIGHT, DOWN, UP };

template<class T>
class geometry_2d_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    // Переменные на основе которых строятся функции форм.
    static constexpr symdiff::variable<0> xi{};
    static constexpr symdiff::variable<1> eta{};

public:
    virtual ~geometry_2d_base() noexcept = default;
    virtual T boundary(const side_2d bound, const T x) const = 0;
};

// Данная реализация подразумевает, что данные о геометрии области наследуются из Shape_Type<T>.
// Таким образом пользователь сможет добавлять свои области не прибегая к дублированию интерфейса.
template<class T, template<class> class Shape_Type>
class geometry_2d : public geometry_2d_base<T>,
                    public Shape_Type<T> {
    static_assert(Shape_Type<T>::boundary.size() == 4, "Wrong number of boundaries.");

public:
    ~geometry_2d() override = default;

    T boundary(const side_2d bound, const T x) const override { return Shape_Type<T>::boundary[size_t(bound)](x); }
};

// Описание классов стратегий Shape_Type<T>.
// Каждый класс описывающий двумерную геометрию должен содержать в себе статический массив из четырёх функий, который называется boundary.
// Каждая функция должна описывать соответствующий предел интегрирования. Естественно такое описание может быть неоднозначным.
// Под неоднозначностью понимается то, что переменные пределы интегрирования могут быть как снизу-сверху, так и слева-справа.
// Выбор подходящего варианта описания зависит от удобства пользования.

template<class T>
class triangle_element_geometry {
public:
    virtual ~triangle_element_geometry() noexcept = default;

protected:
    explicit triangle_element_geometry() = default;
    static inline const std::array<std::function<T(const T)>, 4>
        boundary = { [](const T   ) { return 0.0;    },
                     [](const T   ) { return 1.0;    },
                     [](const T   ) { return 0.0;    },
                     [](const T xi) { return 1.0-xi; } };
};

template<class T>
class rectangle_element_geometry {
public:
    virtual ~rectangle_element_geometry() noexcept = default;

protected:
    explicit rectangle_element_geometry() = default;
    static inline const std::array<std::function<T(const T)>, 4>
        boundary = { [](const T) { return -1.0; },
                     [](const T) { return  1.0; },
                     [](const T) { return -1.0; },
                     [](const T) { return  1.0; } };
};

}

#endif