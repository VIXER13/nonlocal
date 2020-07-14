#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

// В данном модуле описаны различные возмоные функции виляния

#include <cmath>
#include <array>
#include "metamath/power.hpp"

namespace nonlocal::influence {

// Семейства функций с эллиптической областью влияния.
class _elliptical_region {
    _elliptical_region() noexcept = default;

    // Функция вычисления коэффициентов уравнения эллипса для различных радиусов и углов поворота.
    // Угол поворота theta берётся против часовой стрелки.
    template<class Type>
    static std::array<Type, 3> calc_coeff(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        const Type cos_theta = cos(theta), sin_theta = sin(theta);
        return { metamath::power<2>(cos_theta / r[0]) + metamath::power<2>(sin_theta / r[1]),
                                       cos_theta*sin_theta * (2./(r[0]*r[0]) - 2./(r[1]*r[1])),
                 metamath::power<2>(sin_theta / r[0]) + metamath::power<2>(cos_theta / r[1]) };
    }

    // Функция вычисления уравнения эллипса.
    // Если она возвращает значение < 1, то точка находится внутри эллипса, иначе за его пределами.
    template<class Type>
    static Type calc_ellipse(const std::array<Type, 3>& coeff,
                             const std::array<Type, 2>& x, const std::array<Type, 2>& y) noexcept {
        const Type dx = x[0] - x[1], dy = y[0] - y[1];
        return coeff[0]*dx*dx + (coeff[1]*dx + coeff[2]*dy)*dy;
    }

public:
    // Широкий класс полиномиальных функций влияния у которых есть два параметра:
    // параметр p отвечает за равномерность распределения, чем он больше, тем равномернее;
    // параметр q концентрирует влияние в центре области, чем он больше, тем более концентрированное влияние в центре.
    // Оба параметра должны быть больше 0. В наиболее общем случае они так же могут принимать и дробные значения,
    // но было решено отдать предпочтение целочисленным, так как использование дробных параметров
    // во-первых не представляет особого практического интереса, а во-вторых повышает сложность вычислений.
    template<class Type, uint64_t p, uint64_t q>
    friend class polynomial;

    // Константа является предельным случаем полиномиальных функций влияния, в случае, когда p->inf.
    // Была выделена в отдельный класс для ясности намерений и простоты вычислений.
    template<class Type>
    friend class constant;

    // Функция нормального распределния.
    template<class Type>
    friend class normal_distribution;
};

template<class Type, uint64_t p, uint64_t q>
class polynomial {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    static_assert(p > 0, "Parameter p must be greater than 0.");
    static_assert(q > 0, "Parameter q must be greater than 0.");

    std::array<Type, 3> _coeff = {};
    Type _norm = 0;

public:
    polynomial(const Type r) noexcept {
        set_parameters(r);
    }

    polynomial(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        set_parameters(r, theta);
    }

    void set_parameters(const Type r) noexcept {
        set_parameters({r, r});
    }

    void set_parameters(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        _coeff = _elliptical_region::calc_coeff(r, theta);
        _norm = p / ((M_PI + M_PI) * r[0] * r[1] * std::beta(2./p, q+1.));
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        const Type h = _elliptical_region::calc_ellipse(_coeff, {xL, xNL}, {yL, yNL});
        if constexpr (p % 2)
            return h > 1 ? 0. : _norm * metamath::power<q>(1. - metamath::power<p/2>(h)*sqrt(h));
        else
            return h > 1 ? 0. : _norm * metamath::power<q>(1. - metamath::power<p/2>(h));
    }
};

template<class Type>
class constant {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

    std::array<Type, 3> _coeff = {};
    Type _norm = 0;

public:
    constant(const Type r) noexcept {
        set_parameters(r);
    }

    constant(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        set_parameters(r, theta);
    }

    void set_parameters(const Type r) noexcept {
        set_parameters({r, r});
    }

    void set_parameters(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        _coeff = _elliptical_region::calc_coeff(r, theta);
        _norm = 1. / (M_PI * r[0] * r[1]);
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        return _elliptical_region::calc_ellipse(_coeff, {xL, xNL}, {yL, yNL}) < 1. ? _norm : 0.;
    }
};

template<class Type>
class normal_distribution {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

    std::array<Type, 3> _coeff = {};
    Type _norm = 0;

public:
    normal_distribution(const Type r) noexcept {
        set_parameters(r);
    }

    normal_distribution(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        set_parameters(r, theta);
    }

    void set_parameters(const Type r) noexcept {
        set_parameters({r, r});
    }

    void set_parameters(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        _coeff = _elliptical_region::calc_coeff(r, theta);
        _norm = 1. / (M_PI * r[0] * r[1]);
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        return _norm * exp(-_elliptical_region::calc_ellipse(_coeff, {xL, xNL}, {yL, yNL}));
    }
};

}

#endif