#ifndef FINITE_ELEMENT_GEOMETRY_1D_HPP
#define FINITE_ELEMENT_GEOMETRY_1D_HPP

#include <array>

namespace finite_element {

// Одномерную геометрию можно описать началом и концом отрезка.
enum class side_1d : uint8_t {LEFT, RIGHT};

template<class Type>
class geometry_1d_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual Type boundary(const side_1d bound) const = 0;
    virtual ~geometry_1d_base() noexcept = default;
};

// Данная реализация подразумевает, что данные о начале и конце отрезка наследуются из Shape_Type<Type>.
// Таким образом пользователь сможет добавлять свои области не прибегая к дублированию интерфейса.
template<class Type, template<class> class Shape_Type>
class geometry_1d : public geometry_1d_base<Type>,
                    public Shape_Type<Type> {
    static_assert(Shape_Type<Type>::boundary.size() == 2, "Wrong number of boundaries.");

public:
    Type boundary(const side_1d bound) const override { return Shape_Type<Type>::boundary[static_cast<size_t>(bound)]; }
};

// Описание классов стратегий Shape_Type<Type>.
// Каждый класс описывающий одномерную геометрию должен содержать в себе статический массив из двух чисел, который называется boundary.

template<class Type>
class standart_segment_geometry {
public:
    virtual ~standart_segment_geometry() noexcept = default;

protected:
    explicit standart_segment_geometry() noexcept = default;
    static constexpr std::array<Type, 2> boundary = { -1.0, 1.0 };
};

}

#endif