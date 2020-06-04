#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

#include <cmath>
#include <array>
#include "power.hpp"

namespace nonlocal::influence
{

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
        const Type cos_theta = cos(theta), sin_theta = sin(theta);
        _coeff = { math_meta::power<2>(cos_theta / r[0]) + math_meta::power<2>(sin_theta / r[1]),
                                      2.*cos_theta*sin_theta * (1./(r[0]*r[0]) - 1./(r[1]*r[1])),
                   math_meta::power<2>(sin_theta / r[0]) + math_meta::power<2>(cos_theta / r[1]) };
        _norm = 1. / (M_PI * r[0] * r[1]);
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        const Type dx = xL - xNL, dy = yL - yNL;
        return _coeff[0]*dx*dx + (_coeff[1]*dx + _coeff[2]*dy)*dy < 1. ? _norm : 0.;
    }
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
        const Type cos_theta = cos(theta), sin_theta = sin(theta);
        _coeff = { math_meta::power<2>(cos_theta / r[0]) + math_meta::power<2>(sin_theta / r[1]),
                                      2.*cos_theta*sin_theta * (1./(r[0]*r[0]) - 1./(r[1]*r[1])),
                   math_meta::power<2>(sin_theta / r[0]) + math_meta::power<2>(cos_theta / r[1]) };
        _norm = p / ((M_PI + M_PI) * r[0] * r[1] * std::beta(2./p, q+1.));
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        const Type dx = xL - xNL, dy = yL - yNL, h = _coeff[0]*dx*dx + (_coeff[1]*dx + _coeff[2]*dy)*dy;
        if constexpr (p % 2)
            return h > 1 ? 0. : _norm * math_meta::power<q>(1. - math_meta::power<p/2>(h)*sqrt(h));
        else
            return h > 1 ? 0. : _norm * math_meta::power<q>(1. - math_meta::power<p/2>(h));
    }
};

template<class Type>
class normal_distribution {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

    Type _rr = 0, _norm = 0;

public:
    normal_distribution(const Type r) noexcept {
        set_parameters(r);
    }

    void set_parameters(const Type r) noexcept {
        _rr = r*r;
        _norm = 1. / (M_PI * _rr);
    }

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) const noexcept {
        return _norm * exp(-(math_meta::power<2>(xL - xNL) + math_meta::power<2>(yL - yNL)) / _rr);
    }
};

}

#endif