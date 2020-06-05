#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

#include <cmath>
#include <array>
#include "power.hpp"

namespace nonlocal::influence {

class _elliptical_region {
    _elliptical_region() noexcept = default;

    template<class Type>
    static std::array<Type, 3> calc_coeff(const std::array<Type, 2>& r, const Type theta = 0) noexcept {
        const Type cos_theta = cos(theta), sin_theta = sin(theta);
        return { math_meta::power<2>(cos_theta / r[0]) + math_meta::power<2>(sin_theta / r[1]),
                                    2.*cos_theta*sin_theta * (1./(r[0]*r[0]) - 1./(r[1]*r[1])),
                 math_meta::power<2>(sin_theta / r[0]) + math_meta::power<2>(cos_theta / r[1]) };
    }

    template<class Type>
    static Type calc_ellipse(const std::array<Type, 3>& coeff,
                             const std::array<Type, 2>& x, const std::array<Type, 2>& y) noexcept {
        const Type dx = x[0] - x[1], dy = y[0] - y[1];
        return coeff[0]*dx*dx + (coeff[1]*dx + coeff[2]*dy)*dy;
    }

public:
    template<class Type>
    friend class constant;

    template<class Type, uint64_t p, uint64_t q>
    friend class polynomial;

    template<class Type>
    friend class normal_distribution;
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
            return h > 1 ? 0. : _norm * math_meta::power<q>(1. - math_meta::power<p/2>(h)*sqrt(h));
        else
            return h > 1 ? 0. : _norm * math_meta::power<q>(1. - math_meta::power<p/2>(h));
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