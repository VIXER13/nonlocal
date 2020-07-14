#ifndef FINITE_ELEMENT_GEOMETRY_1D_HPP
#define FINITE_ELEMENT_GEOMETRY_1D_HPP

#include <array>
#include "../../symdiff/symdiff.hpp"

namespace metamath::finite_element {

// Одномерную геометрию можно описать началом и концом отрезка.
enum class side_1d : uint8_t { LEFT, RIGHT };

template<class T>
class geometry_1d_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    // Переменная, на основе которой строятся функции формы элементов.
    static constexpr symdiff::variable<0> xi{};

public:
    virtual ~geometry_1d_base() noexcept = default;
    virtual T boundary(const side_1d bound) const = 0;
};

// Данная реализация подразумевает, что данные о начале и конце отрезка наследуются из Shape_Type<T>.
// Таким образом пользователь сможет добавлять свои области не прибегая к дублированию интерфейса.
template<class T, template<class> class Shape_Type>
class geometry_1d : public geometry_1d_base<T>,
                    public Shape_Type<T> {
    static_assert(Shape_Type<T>::boundary.size() == 2, "Wrong number of boundaries.");

public:
    T boundary(const side_1d bound) const override { return Shape_Type<T>::boundary[size_t(bound)]; }
};

// Описание классов стратегий Shape_Type<T>.
// Каждый класс описывающий одномерную геометрию должен содержать в себе статический массив из двух чисел, который называется boundary.

template<class T>
class standart_segment_geometry {
public:
    virtual ~standart_segment_geometry() noexcept = default;

protected:
    explicit standart_segment_geometry() noexcept = default;
    static constexpr std::array<T, 2> boundary = { -1.0, 1.0 };
};

}

#endif