#pragma once

#include <metamath/metamath.hpp>
#include "../../equation_parameters.hpp"

#include <string>
#include <unordered_map>
#include <variant>

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

    using _base = parameter_2d_base<T>;

public:
    using conductivity_t = 
        std::conditional_t<
            Coefficients == coefficients_t::SOLUTION_DEPENDENT,
            std::function<T(const std::array<T, 2>&, const T)>,
            std::conditional_t<
                Coefficients == coefficients_t::SPACE_DEPENDENT,
                std::function<T(const std::array<T, 2>&)>, T
            >
        >;

    using conductivity_variants = std::variant<
        conductivity_t, 
        std::array<conductivity_t, 2>,
        metamath::types::square_matrix<conductivity_t, 2>
    >;

private:
    static constexpr metamath::types::square_matrix<conductivity_t, 2> init() noexcept {
        if constexpr (Coefficients == coefficients_t::CONSTANTS)
            return {T{1}, T{0}, T{0}, T{0}};
        if constexpr (Coefficients == coefficients_t::SPACE_DEPENDENT)
            return {
                [](const std::array<T, 2>& x) constexpr noexcept { return T{1}; },
                [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; },
                [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; },
                [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
            };
        if constexpr (Coefficients == coefficients_t::SOLUTION_DEPENDENT)
            return {
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return T{1}; },
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return T{0}; },
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return T{0}; },
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return T{0}; }
            };
    }

    material_t get_coefficient(const conductivity_variants& conductivity) {
        if (std::holds_alternative<conductivity_t>(conductivity))
            return material_t::ISOTROPIC;
        if (std::holds_alternative<std::array<conductivity_t, 2>>(conductivity))
            return material_t::ORTHOTROPIC;
        return material_t::ANISOTROPIC;
    }

public:
    metamath::types::square_matrix<conductivity_t, 2> conductivity = init();

    constexpr parameter_2d() noexcept
        : parameter_2d_base<T>{Coefficients} {}
    explicit parameter_2d(const conductivity_variants& cond, const T capacity = T{1},
                          const T density = T{1}) noexcept(Coefficients == coefficients_t::CONSTANTS)
        : parameter_2d_base<T>{Coefficients, capacity, density, get_coefficient(cond)} {
        if (_base::material == material_t::ISOTROPIC)
            conductivity[X][X] = std::get<conductivity_t>(cond);
        else if (_base::material == material_t::ORTHOTROPIC) {
            const auto& tmp = std::get<std::array<conductivity_t, 2>>(cond);
            conductivity[X][X] = tmp[X];
            conductivity[Y][Y] = tmp[Y];
        } else
            conductivity = std::get<metamath::types::square_matrix<conductivity_t, 2>>(cond);
    }
    ~parameter_2d() noexcept override = default;
};

template<class T>
using parameter_2d_sptr = std::shared_ptr<parameter_2d_base<T>>;

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d_sptr>>;

template<coefficients_t Coefficients, bool Check_Nullptr = true, class T>
constexpr parameter_2d<T, Coefficients>* parameter_cast(parameter_2d_base<T> *const parameter) noexcept(Check_Nullptr) {
    if constexpr (Check_Nullptr)
        return parameter && parameter->type == Coefficients ? 
            static_cast<parameter_2d<T, Coefficients>*>(parameter) : nullptr;
    else
        return parameter->type == Coefficients ? 
            static_cast<parameter_2d<T, Coefficients>*>(parameter) : nullptr;
}

template<coefficients_t Coefficients, bool Check_Nullptr = true, class T>
constexpr const parameter_2d<T, Coefficients>* parameter_cast(const parameter_2d_base<T> *const parameter) noexcept(Check_Nullptr) {
    if constexpr (Check_Nullptr)
        return parameter && parameter->type == Coefficients ? 
            static_cast<const parameter_2d<T, Coefficients>*>(parameter) : nullptr;
    else
        return parameter->type == Coefficients ? 
            static_cast<parameter_2d<T, Coefficients>*>(parameter) : nullptr;
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

};