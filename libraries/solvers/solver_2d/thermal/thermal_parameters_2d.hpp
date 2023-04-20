#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

#include "metamath.hpp"
#include "../../equation_parameters.hpp"

#include <string>
#include <unordered_map>

namespace nonlocal::thermal {

template<class T>
struct parameter_2d_base {
    const coefficients_t type;

    T capacity = T{1};
    T density = T{1};
    material_t material = material_t::ISOTROPIC;

protected:
    explicit parameter_2d_base(const coefficients_t coefficients, const T capacity = T{1}, 
                               const T density = T{1}, const material_t material = material_t::ISOTROPIC) noexcept
        : type{coefficients}
        , capacity{capacity}
        , density{density}
        , material{material} {}

public:
    virtual ~parameter_2d_base() noexcept = default;
};

template<class T, coefficients_t Coefficients = coefficients_t::CONSTANTS>
class parameter_2d final : public parameter_2d_base<T> {
    static_assert(Coefficients == coefficients_t::CONSTANTS ||
                  Coefficients == coefficients_t::SPACE_DEPENDENT ||
                  Coefficients == coefficients_t::SOLUTION_DEPENDENT,
                  "Unknown coefficient type.");

public:
    using conductivity_t = metamath::types::square_matrix<
        std::conditional_t<
            Coefficients == coefficients_t::SOLUTION_DEPENDENT,
            std::function<T(const T, const T)>,
            std::conditional_t<
                Coefficients == coefficients_t::SPACE_DEPENDENT,
                std::function<T(const T)>, T
            >
        >,
        2
    >;

private:
    static constexpr conductivity_t init() noexcept {
        if constexpr (Coefficients == coefficients_t::CONSTANTS)
            return {T{1}};
        if constexpr (Coefficients == coefficients_t::SPACE_DEPENDENT)
            return {
                [](const T x) constexpr noexcept { return T{1}; },
                [](const T x) constexpr noexcept { return T{0}; },
                [](const T x) constexpr noexcept { return T{0}; },
                [](const T x) constexpr noexcept { return T{0}; }
            };
        if constexpr (Coefficients == coefficients_t::SOLUTION_DEPENDENT)
            return {
                [](const T x, const T solution) constexpr noexcept { return T{1}; },
                [](const T x, const T solution) constexpr noexcept { return T{0}; },
                [](const T x, const T solution) constexpr noexcept { return T{0}; },
                [](const T x, const T solution) constexpr noexcept { return T{0}; }
            };
    }

public:
    conductivity_t conductivity = init();

    constexpr parameter_2d() noexcept
        : parameter_2d_base<T>{Coefficients} {}
    explicit parameter_2d(const conductivity_t& conductivity, const T capacity = T{1}, 
                          const T density = T{1}, const material_t material = material_t::ISOTROPIC) noexcept(Coefficients == coefficients_t::CONSTANTS)
        : parameter_2d_base<T>{Coefficients, capacity, density, material}
        , conductivity{conductivity} {}
    ~parameter_2d() noexcept override = default;
};

template<class T>
using parameter_2d_sptr = std::shared_ptr<parameter_2d_base<T>>;

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d_sptr>>;

template<coefficients_t Coefficients, class T>
constexpr parameter_2d<T, Coefficients>* parameter_cast(parameter_2d_base<T> *const parameter) noexcept {
    return parameter && parameter->type == Coefficients ? 
        static_cast<parameter_2d<T, Coefficients>*>(parameter) : nullptr;
}

template<coefficients_t Coefficients, class T>
constexpr const parameter_2d<T, Coefficients>* parameter_cast(const parameter_2d_base<T> *const parameter) noexcept {
    return parameter && parameter->type == Coefficients ? 
        static_cast<const parameter_2d<T, Coefficients>*>(parameter) : nullptr;
}

template<coefficients_t Coefficients, class T>
parameter_2d<T, Coefficients>& parameter_cast(parameter_2d_base<T>& parameter) {
    if (parameter.type == Coefficients)
        return static_cast<parameter_2d<T, Coefficients>&>(parameter);
    throw std::domain_error{"Can't cast parameter."};
}

template<coefficients_t Coefficients, class T>
const parameter_2d<T, Coefficients>& parameter_cast(const parameter_2d_base<T>& parameter) {
    if (parameter.type == Coefficients)
        return static_cast<const parameter_2d<T, Coefficients>&>(parameter);
    throw std::domain_error{"Can't cast parameter."};
}

//template<class T>
//struct parameter_2d final {
//    metamath::types::square_matrix<T, 2> conductivity = {T{1}};
//    T capacity = T{1};
//    T density = T{1};
//    material_t material = material_t::ISOTROPIC;
//};

//template<class T>
//using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

};

#endif