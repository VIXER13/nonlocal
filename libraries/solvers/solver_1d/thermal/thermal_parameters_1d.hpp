#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

#include "../../equation_parameters.hpp"
#include <memory>

namespace nonlocal::thermal {

template<class T>
struct parameter_1d_base {
    const coefficients_t type;
    T capacity = T{1};
    T density = T{1};

protected:
    explicit parameter_1d_base(const coefficients_t coefficients, const T capacity = T{1}, const T density = T{1}) noexcept
        : type{coefficients}
        , capacity{capacity}
        , density{density} {}

public:
    virtual ~parameter_1d_base() noexcept = default;
};

template<class T, coefficients_t Coefficients = coefficients_t::CONSTANTS>
class parameter_1d final : public parameter_1d_base<T> {
    static_assert(Coefficients == coefficients_t::CONSTANTS ||
                  Coefficients == coefficients_t::SPACE_DEPENDENT ||
                  Coefficients == coefficients_t::SOLUTION_DEPENDENT,
                  "Unknown coefficient type.");

    static constexpr auto init() noexcept {
        if constexpr (Coefficients == coefficients_t::CONSTANTS)
            return T{1};
        if constexpr (Coefficients == coefficients_t::SPACE_DEPENDENT)
            return [](const T x) constexpr noexcept { return T{1}; };
        if constexpr (Coefficients == coefficients_t::SOLUTION_DEPENDENT)
            return [](const T x, const T solution) constexpr noexcept { return T{1}; };
    }

public:
    using conductivity_t = 
        std::conditional_t<
            Coefficients == coefficients_t::SOLUTION_DEPENDENT,
            std::function<T(const T, const T)>,
            std::conditional_t<
                Coefficients == coefficients_t::SPACE_DEPENDENT,
                std::function<T(const T)>, T
            >
        >;

    conductivity_t conductivity = init();

    constexpr parameter_1d() noexcept
        : parameter_1d_base<T>{Coefficients} {}
    explicit parameter_1d(const conductivity_t& conductivity, const T capacity = T{1}, const T density = T{1}) noexcept(Coefficients == coefficients_t::CONSTANTS)
        : parameter_1d_base<T>{Coefficients, capacity, density}
        , conductivity{conductivity} {}
    ~parameter_1d() noexcept override = default;
};

template<class T>
using parameter_1d_sptr = std::shared_ptr<parameter_1d_base<T>>;

template<class T>
using parameters_1d = std::vector<equation_parameters<1, T, parameter_1d_sptr>>;

template<coefficients_t Coefficients, class T>
constexpr parameter_1d<T, Coefficients>* parameter_cast(parameter_1d_base<T> *const parameter) noexcept {
    return parameter && parameter->type == Coefficients ? 
        static_cast<parameter_1d<T, Coefficients>*>(parameter) : nullptr;
}

template<coefficients_t Coefficients, class T>
constexpr const parameter_1d<T, Coefficients>* parameter_cast(const parameter_1d_base<T> *const parameter) noexcept {
    return parameter && parameter->type == Coefficients ? 
        static_cast<const parameter_1d<T, Coefficients>*>(parameter) : nullptr;
}

template<coefficients_t Coefficients, class T>
parameter_1d<T, Coefficients>& parameter_cast(parameter_1d_base<T>& parameter) {
    if (parameter.type == Coefficients)
        return static_cast<parameter_1d<T, Coefficients>&>(parameter);
    throw std::domain_error{"Can't cast parameter."};
}

template<coefficients_t Coefficients, class T>
const parameter_1d<T, Coefficients>& parameter_cast(const parameter_1d_base<T>& parameter) {
    if (parameter.type == Coefficients)
        return static_cast<const parameter_1d<T, Coefficients>&>(parameter);
    throw std::domain_error{"Can't cast parameter."};
}

}

#endif