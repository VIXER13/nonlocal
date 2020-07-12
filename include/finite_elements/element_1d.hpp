#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

// Классы стратегий для одномерных элементов.

#include "geometry_1d.hpp"
#include <functional>

namespace finite_element {

// Ниже представлены некоторые реализации классов ElementType.
// Данный список можно пополнить.

// Ниже описаны элементы с базисами в виде полиномов Лагранжа.

template<class Type>
class linear : protected geometry_1d<Type, standart_segment_geometry> {
protected:
    explicit linear() noexcept = default;

    // Нумерация узлов на линейном элементе: 0--1
    static inline constexpr std::array<Type, 2> nodes = { -1.0, 1.0 };
    static inline const std::array<std::function<Type(const Type)>, 2>
        basic_N   = { [](const Type xi) { return 0.5*(1.0-xi); },
                      [](const Type xi) { return 0.5*(1.0+xi); } },

        basic_Nxi = { [](const Type   ) { return -0.5; },
                      [](const Type   ) { return  0.5; } };
};

template<class Type>
class quadratic : protected geometry_1d<Type, standart_segment_geometry> {
protected:
    explicit quadratic() noexcept = default;

    // Нумерация узлов на квадратичном элементе: 0--1--2
    static inline constexpr std::array<Type, 3> nodes = { -1.0, 0.0, 1.0 };
    static inline const std::array<std::function<Type(const Type)>, 3>
        basic_N   = { [](const Type xi) { return  -0.5*xi*(1.0-xi); },
                      [](const Type xi) { return (1.0+xi)*(1.0-xi); },
                      [](const Type xi) { return   0.5*xi*(1.0+xi); } },

        basic_Nxi = { [](const Type xi) { return -0.5+xi; },
                      [](const Type xi) { return -2.0*xi; },
                      [](const Type xi) { return  0.5+xi; } };;
};

template<class Type>
class qubic : protected geometry_1d<Type, standart_segment_geometry> {
protected:
    explicit qubic() noexcept = default;

    // Нумерация узлов на кубическом элементе: 0--1--2--3
    static inline constexpr std::array<Type, 4> nodes = { -1.0, -1.0/3.0, 1.0/3.0, 1.0 };
    static inline const std::array<std::function<Type(const Type)>, 4>
        basic_N   = { [](const Type xi) { return -0.5625*(xi-1.0)*(xi+1.0/3.0)*(xi-1.0/3.0); },
                      [](const Type xi) { return  1.6875*(xi+1.0)*(xi-1.0)    *(xi-1.0/3.0); },
                      [](const Type xi) { return -1.6875*(xi+1.0)*(xi-1.0)    *(xi+1.0/3.0); },
                      [](const Type xi) { return  0.5625*(xi+1.0)*(xi+1.0/3.0)*(xi-1.0/3.0); } },
                     
        basic_Nxi = { [](const Type xi) { return (-1.6875*xi + 1.125)*xi + 0.0625; },
                      [](const Type xi) { return ( 5.0625*xi - 1.125)*xi - 1.6875; },
                      [](const Type xi) { return (-5.0625*xi - 1.125)*xi + 1.6875; },
                      [](const Type xi) { return ( 1.6875*xi + 1.125)*xi - 0.0625; } };
};

}

#endif