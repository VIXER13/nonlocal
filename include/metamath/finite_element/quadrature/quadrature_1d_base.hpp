#ifndef FINITE_ELEMENT_QUADRATURE_BASE_1D_HPP
#define FINITE_ELEMENT_QUADRATURE_BASE_1D_HPP

#include "quadrature_base.hpp"
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class T>
class quadrature_1d_base : public quadrature_base<T> {
public:
    ~quadrature_1d_base() override = default;
    virtual const std::array<T, 1>& node(const size_t i) const = 0; // Получение координаты узла под номером i
    virtual T boundary(const side_1d bound) const = 0; // Геометрия квадратуры.
};

}

#endif