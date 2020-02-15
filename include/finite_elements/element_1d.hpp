#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

// Реализация класса одномерных элементов.

#include "element_base.hpp"

namespace finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса ElementType. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class Type, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<Type>,
                   public Element_Type<Type> {
    static_assert(Element_Type<Type>::basicN.size() == Element_Type<Type>::basicNxi.size(),
                  "The number of functions and their derivatives does not match.");
public:
    size_t nodes_count() const override { return Element_Type<Type>::basicN.size(); }
        
    Type N  (const size_t i, const Type xi) const override { return Element_Type<Type>::basicN  [i](xi); }
    Type Nxi(const size_t i, const Type xi) const override { return Element_Type<Type>::basicNxi[i](xi); }

    Type boundary(const side_1d bound) const override { return Element_Type<Type>::boundary(bound); }

    virtual ~element_1d() = default;
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Ниже представлены некоторые реализации классов ElementType.
// Данный список можно пополнить.

// Ниже описаны элементы с базисами в виде полиномов Лагранжа.

template<class Type>
class linear : protected geometry_1d<Type, standart_segment_geometry> {
public:
    virtual ~linear() = default;

protected:
    explicit linear() noexcept = default;
    // Нумерация узлов на линейном элементе: 0--1
    static inline const std::array<std::function<Type(const Type)>, 2>
        basicN   = { [](const Type xi) { return 0.5*(1.0-xi); },
                     [](const Type xi) { return 0.5*(1.0+xi); } },

        basicNxi = { [](const Type xi) { return -0.5; },
                     [](const Type xi) { return  0.5; } };
};

template<class Type>
class quadratic : protected geometry_1d<Type, standart_segment_geometry> {
public:
    virtual ~quadratic() = default;

protected:
    explicit quadratic() noexcept = default;
    // Нумерация узлов на квадратичном элементе: 0--1--2
    static inline const std::array<std::function<Type(const Type)>, 3>
        basicN   = { [](const Type xi) { return  -0.5*xi*(1.0-xi); },
                     [](const Type xi) { return (1.0+xi)*(1.0-xi); },
                     [](const Type xi) { return   0.5*xi*(1.0+xi); } },

        basicNxi = { [](const Type xi) { return -0.5+xi; },
                     [](const Type xi) { return -2.0*xi; },
                     [](const Type xi) { return  0.5+xi; } };;
};

template<class Type>
class qubic : protected geometry_1d<Type, standart_segment_geometry> {
public:
    virtual ~qubic() = default;
    
protected:
    explicit qubic() noexcept = default;
    // Нумерация узлов на кубическом элементе: 0--1--2--3
    static inline const std::array<std::function<Type(const Type)>, 4>
        basicN   = { [](const Type xi) { return -0.5625*(xi-1.0)*(xi+1.0/3.0)*(xi-1.0/3.0); },
                     [](const Type xi) { return  1.6875*(xi+1.0)*(xi-1.0)    *(xi-1.0/3.0); },
                     [](const Type xi) { return -1.6875*(xi+1.0)*(xi-1.0)    *(xi+1.0/3.0); },
                     [](const Type xi) { return  0.5625*(xi+1.0)*(xi+1.0/3.0)*(xi-1.0/3.0); } },
                     
        basicNxi = { [](const Type xi) { return (-1.6875*xi + 1.125)*xi + 0.0625; },
                     [](const Type xi) { return ( 5.0625*xi - 1.125)*xi - 1.6875; },
                     [](const Type xi) { return (-5.0625*xi - 1.125)*xi + 1.6875; },
                     [](const Type xi) { return ( 1.6875*xi + 1.125)*xi - 0.0625; } };
};

#pragma GCC diagnostic pop

}

#endif