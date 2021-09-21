#ifndef FINITE_ELEMENT_QUADRATURE_1D_HPP
#define FINITE_ELEMENT_QUADRATURE_1D_HPP

#include "quadrature_1d_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные об узлах и весах квадратуры, а так же геометрия наследуются от класса Quadrature_Type.
// Таким образом пользователь сможет добавлять свои квадратуры не прибегая к дублированию интерфейса.
template<class T, template<class> class Quadrature_Type>
class quadrature_1d : public quadrature_1d_base<T>,
                      public Quadrature_Type<T> {
    static_assert(Quadrature_Type<T>::nodes.size() == Quadrature_Type<T>::weights.size(),
                  "The number of nodes and weights does not match.");

public:
    ~quadrature_1d() override = default;

    std::unique_ptr<quadrature_1d_base<T>> clone() const override { return std::make_unique<quadrature_1d<T, Quadrature_Type>>(); }

    size_t nodes_count() const override { return Quadrature_Type<T>::nodes.size(); }

    const std::array<T, 1>& node (const size_t i) const override { return Quadrature_Type<T>::nodes[i]; }
    T weight(const size_t i) const override { return Quadrature_Type<T>::weights[i]; }

    T boundary(const side_1d bound) const override { return Quadrature_Type<T>::boundary(bound); }
};

}

#endif