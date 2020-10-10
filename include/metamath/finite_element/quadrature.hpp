#ifndef FINITE_ELEMENT_QUADRATURE_HPP
#define FINITE_ELEMENT_QUADRATURE_HPP

// В данном модуле описан базовый интерфейс и реализация класса квадратуры,
// а так же небольшая коллекция гауссовых квадратур.

#include <cmath>
#include "element_1d_strategy/geometry_1d.hpp"

namespace metamath::finite_element {

template<class T>
class quadrature_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    virtual size_t nodes_count() const = 0; // Количество узлов квадратуры

    virtual T node(const size_t i) const = 0; // Получение координаты узла под номером i
    virtual T weight(const size_t i) const = 0; // Получение веса для узла под номером i

    virtual T boundary(const side_1d bound) const = 0; // Геометрия квадратуры.

    virtual ~quadrature_base() noexcept = default;
};

// Данная реализация подразумевает, что данные об узлах и весах квадратуры, а так же геометрия наследуются от класса Quadrature_Type.
// Таким образом пользователь сможет добавлять свои квадратуры не прибегая к дублированию интерфейса.
template<class T, template<class> class Quadrature_Type>
class quadrature : public quadrature_base<T>,
                   public Quadrature_Type<T> {
    static_assert(Quadrature_Type<T>::nodes.size() == Quadrature_Type<T>::weights.size(),
                  "The number of nodes and weights does not match.");

public:
    ~quadrature() override = default;

    size_t nodes_count() const override { return Quadrature_Type<T>::nodes.size(); }
        
    T node  (const size_t i) const override { return Quadrature_Type<T>::nodes[i]; }
    T weight(const size_t i) const override { return Quadrature_Type<T>::weights[i]; }

    T boundary(const side_1d bound) const override { return Quadrature_Type<T>::boundary(bound); }
};

// Ниже представлены некоторые реализации классов Quadrature_Type.
// Данный список можно пополнить.

// Наследование квадратур от класса геометрии подразумевает возможность использования нестандартных квадратур,
// а так же многомерных квадратур, которые не получаются путём декартова произведения одномерных квадратур.

template<class T>
class gauss1 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss1() noexcept = default;
    ~gauss1() override = default;

    static constexpr std::array<T, 1> nodes   = { 0. },
                                      weights = { 2. };
};

template<class T>
class gauss2 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss2() noexcept = default;
    ~gauss2() override = default;

    static constexpr std::array<T, 2> nodes   = { -1./sqrt(3.), 1./sqrt(3.) }, 
                                      weights = {           1.,          1. };
};

template<class T>
class gauss3 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss3() noexcept = default;
    ~gauss3() override = default;

    static constexpr std::array<T, 3> nodes   = { -sqrt(0.6),    0., sqrt(0.6) }, 
                                      weights = {      5./9., 8./9.,     5./9. };
};

template<class T>
class gauss4 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss4() noexcept = default;
    ~gauss4() override = default;

    static constexpr std::array<T, 4> nodes   = { -sqrt(3./7. + 2./7.*sqrt(1.2)),
                                                  -sqrt(3./7. - 2./7.*sqrt(1.2)),
                                                   sqrt(3./7. - 2./7.*sqrt(1.2)),
                                                   sqrt(3./7. + 2./7.*sqrt(1.2)) }, 

                                      weights = { (18. - sqrt(30.)) / 36.,
                                                  (18. + sqrt(30.)) / 36.,
                                                  (18. + sqrt(30.)) / 36.,
                                                  (18. - sqrt(30.)) / 36. };
};

template<class T>
class gauss5 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss5() noexcept = default;
    ~gauss5() override = default;
    
    static constexpr std::array<T, 5> nodes   = { -1./3. * sqrt(5. + 2.*sqrt(10./7.)),
                                                  -1./3. * sqrt(5. - 2.*sqrt(10./7.)),
                                                                                   0.,
                                                   1./3. * sqrt(5. - 2.*sqrt(10./7.)),
                                                   1./3. * sqrt(5. + 2.*sqrt(10./7.)) }, 

                                      weights = { (322. - 13.*sqrt(70.)) / 900.,
                                                  (322. + 13.*sqrt(70.)) / 900.,
                                                                    128. / 225.,
                                                  (322. + 13.*sqrt(70.)) / 900.,
                                                  (322. - 13.*sqrt(70.)) / 900. };
};

}

#endif