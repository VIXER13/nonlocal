#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

#include <cmath>
#include "power.hpp"

namespace influence_function
{

template<class Type>
class constant
{
    static_assert(std::is_floating_point<Type>::value, "Type must be floating point.");

    const Type r, normalization;

public:
    constant(const Type r) :
        r(r), normalization(1. / (M_PI * r * r)) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL)
    {
        const Type dx = xL-xNL, dy = yL-yNL;
        return sqrt(dx*dx + dy*dy) < r ? normalization : 0.;
    }
};

template<class Type, int64_t p, int64_t q>
class polinomial
{
    static_assert(std::is_floating_point<Type>::value, "Type must be floating point.");
    static_assert(p > 0, "Parameter p must be greater than 0.");
    static_assert(q > 0, "Parameter q must be greater than 0.");

    const Type r, normalization;

public:
    polinomial(const Type r) :
        r(r), normalization(p / (2.*M_PI * r*r * std::beta(2./p, q+1.))) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL)
    {
        const Type dx = xL-xNL, dy = yL-yNL, h = sqrt(dx*dx + dy*dy);
        return h < r ? normalization * math_meta::power<q>(1. - math_meta::power<p>(h/r)) : 0.;
    }
};

template<class Type>
class normal_distribution
{
    static_assert(std::is_floating_point<Type>::value, "Type must be floating point.");

    const Type r_square, normalization;

public:
    normal_distribution(const Type r) :
        r_square(r*r), normalization(1. / (M_PI * r * r)) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL)
    {
        const Type dx = xL-xNL, dy = yL-yNL;
        return normalization * exp(-(dx*dx + dy*dy) / r_square);
    }
};

}

#endif