#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

#include <cmath>
#include "power.hpp"

namespace influence_function
{

template<class Type>
class constant {
    static_assert(std::is_floating_point<Type>::value, "The Type must be floating point.");

    const Type r, normalization;

public:
    constant(const Type r) :
        r(r), normalization(1. / (M_PI * r * r)) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) {
        return sqrt(math_meta::power<2>(xL-xNL) + math_meta::power<2>(yL-yNL)) < r ? normalization : 0.;
    }
};

template<class Type, uint64_t p, uint64_t q>
class polinomial {
    static_assert(std::is_floating_point<Type>::value, "The Type must be floating point.");
    static_assert(p > 0, "Parameter p must be greater than 0.");
    static_assert(q > 0, "Parameter q must be greater than 0.");

    const Type r, normalization;

public:
    polinomial(const Type r) :
        r(r), normalization(p / (2.*M_PI * r*r * std::beta(2./p, q+1.))) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) {
        const Type h = sqrt(math_meta::power<2>(xL-xNL) + math_meta::power<2>(yL-yNL));
        return h < r ? normalization * math_meta::power<q>(1. - math_meta::power<p>(h/r)) : 0.;
    }
};

template<class Type>
class normal_distribution {
    static_assert(std::is_floating_point<Type>::value, "The Type must be floating point.");

    const Type r_square, normalization;

public:
    normal_distribution(const Type r) :
        r_square(r*r), normalization(1. / (M_PI * r * r)) {}

    Type operator()(Type xL, Type xNL, Type yL, Type yNL) {
        return normalization * exp(-(math_meta::power<2>(xL-xNL) + math_meta::power<2>(yL-yNL)) / r_square);
    }
};

}

#endif