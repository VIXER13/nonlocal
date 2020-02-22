#ifndef FINITE_ELEMENT_GEOMETRY_2D_HPP
#define FINITE_ELEMENT_GEOMETRY_2D_HPP

#include <array>
#include <functional>

namespace finite_element {

// Двумерную геометрию можно описать четырьями функциями, каждая из которых описывает границу интегрирования по каждой из сторон.
enum class side_2d : uint8_t {LEFT, RIGHT, DOWN, UP};

template<class Type>
class geometry_2d_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual Type boundary(const side_2d bound, const Type x) const = 0;

    virtual ~geometry_2d_base() noexcept = default;
};

// Данная реализация подразумевает, что данные о геометрии области наследуются из Shape_Type<Type>.
// Таким образом пользователь сможет добавлять свои области не прибегая к дублированию интерфейса.
template<class Type, template<class> class Shape_Type>
class geometry_2d : public geometry_2d_base<Type>,
                    public Shape_Type<Type> {
    static_assert(Shape_Type<Type>::boundary.size() == 4, "Wrong number of boundaries.");

public:
    Type boundary(const side_2d bound, const Type x) const override { return Shape_Type<Type>::boundary[static_cast<size_t>(bound)](x); }
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Описание классов стратегий Shape_Type<Type>.
// Каждый класс описывающий двумерную геометрию должен содержать в себе статический массив из четырёх функий, который называется boundary.
// Каждая функция должна описывать соответствующий предел интегрирования. Естественно такое описание может быть неоднозначным.
// Под неоднозначностью понимается то, что переменные пределы интегрирования могут быть как снизу-сверху, так и слева-справа.
// Выбор подходящего варианта описания зависит от удобства пользования.

template<class Type>
class triangle_element_geometry {
public:
    virtual ~triangle_element_geometry() noexcept = default;

protected:
    explicit triangle_element_geometry() noexcept = default;
    static inline const std::array<std::function<Type(const Type)>, 4>
        boundary = { [](const Type eta) { return 0.0;    },
                     [](const Type eta) { return 1.0;    },
                     [](const Type  xi) { return 0.0;    },
                     [](const Type  xi) { return 1.0-xi; } };
};

template<class Type>
class rectangle_element_geometry {
public:
    virtual ~rectangle_element_geometry() noexcept = default;

protected:
    explicit rectangle_element_geometry() noexcept = default;
    static inline const std::array<std::function<Type(const Type)>, 4>
        boundary = { [](const Type eta) { return -1.0; },
                     [](const Type eta) { return  1.0; },
                     [](const Type  xi) { return -1.0; },
                     [](const Type  xi) { return  1.0; } };
};

#pragma GCC diagnostic pop

}

#endif