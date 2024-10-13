#ifndef NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP
#define NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP

#include "metamath.hpp"
#include "../../equation_parameters.hpp"

#include <array>
#include <string>
#include <unordered_map>

namespace nonlocal::mechanical {

enum class plane_t : bool {
    STRESS,
    STRAIN
};


template<class T>
using hooke_matrix = std::array<T, 6>;

// Isotropic case
// E / (1 - nu * nu)       E * nu / (1 - nu * nu)  0
// E * nu / (1 - nu * nu)  E / (1 - nu * nu)       0
// 0                       0                       E / (2 * (1 + nu))

// Orthotropic case
// Ex / (1 - nuXY * nuYX)         Ey * nuXY / (1 - nuXY * nuYX)  0
// Ex * nuYX / (1 - nuXY * nuYX)  Ey / (1 - nuXY * nuYX)         0
// 0                              0                              Gxy

template<class T>
struct parameter_2d_base {
    const coefficients_t type;

    material_t material = material_t::ISOTROPIC;

protected:
    explicit parameter_2d_base(const coefficients_t coefficients, const material_t material = material_t::ISOTROPIC) noexcept
        : type{coefficients}
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
    using dependency_t = 
        std::conditional_t<
            Coefficients == coefficients_t::SOLUTION_DEPENDENT,
            std::function<T(const std::array<T, 2>&, const T)>,
            std::conditional_t<
                Coefficients == coefficients_t::SPACE_DEPENDENT,
                std::function<T(const std::array<T, 2>&)>, T
            >
        >;

private:
    static constexpr std::array<dependency_t, 2> init(const std::array<T, 2>& initial_coefficients) noexcept {
        if constexpr (Coefficients == coefficients_t::CONSTANTS)
            return {initial_coefficients[0], initial_coefficients[1]};
        if constexpr (Coefficients == coefficients_t::SPACE_DEPENDENT)
            return {
                [](const std::array<T, 2>& x) constexpr noexcept { return initial_coefficients[0]; },
                [](const std::array<T, 2>& x) constexpr noexcept { return initial_coefficients[1]; }
            };
        if constexpr (Coefficients == coefficients_t::SOLUTION_DEPENDENT)
            return {
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return initial_coefficients[0]; },
                [](const std::array<T, 2>& x, const T solution) constexpr noexcept { return initial_coefficients[1]; }
            };
    }

public:
    std::array<dependency_t, 2> youngs_modulus = init(std::array<T, 2>{T{210}, T{210}});
    std::array<dependency_t, 2> poissons_ratio = init(std::array<T, 2>{T{0.3}, T{0.3}});
    dependency_t shear_modulus_xy = 210;
    T thermal_expansion = 13e-6;
        
    constexpr T E(const plane_t plane) const noexcept;
    constexpr T nu(const plane_t plane) const noexcept;
    constexpr hooke_matrix<T> hooke(const plane_t plane) const noexcept;

    constexpr parameter_2d() noexcept
        : parameter_2d_base<T>{Coefficients} {}
    explicit parameter_2d(const dependency_t& youngs_modulus, const dependency_t& poissons_ratio) 
                                                              noexcept(Coefficients == coefficients_t::CONSTANTS)
        : parameter_2d_base<T>{Coefficients,  material_t::ISOTROPIC}
        , youngs_modulus{youngs_modulus}, poissons_ratio{poissons_ratio} {}
    explicit parameter_2d(const std::array<dependency_t, 2>& youngs_modulus, const std::array<dependency_t, 2>& poissons_ratio,
                         const dependency_t shear_modulus_xy) noexcept(Coefficients == coefficients_t::CONSTANTS)
        : parameter_2d_base<T>{Coefficients, material_t::ORTHOTROPIC}
        , youngs_modulus{youngs_modulus}, poissons_ratio{poissons_ratio}, shear_modulus_xy{shear_modulus_xy}  {}
    ~parameter_2d() noexcept override = default;
};

template<class T>
using parameter_2d_sptr = std::shared_ptr<parameter_2d_base<T>>;

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d_sptr>>;

template<class T>
struct mechanical_parameters_2d final {
    parameters_2d<T> materials;
    std::vector<T> delta_temperature; // if empty, then thermal expansions are not taken into account
    plane_t plane = plane_t::STRESS;
};

template<class T, coefficients_t Coefficients>
constexpr T parameter_2d<T, Coefficients>::E(const plane_t plane) const noexcept {
    return plane == plane_t::STRESS ? youngs_modulus : youngs_modulus / (T{1} - poissons_ratio * poissons_ratio);
}

template<class T, coefficients_t Coefficients>
constexpr T parameter_2d<T, Coefficients>::nu(const plane_t plane) const noexcept {
    return plane == plane_t::STRESS ? poissons_ratio : poissons_ratio / (T{1} - poissons_ratio);
}

template<class T, coefficients_t Coefficients>
constexpr hooke_matrix<T> parameter_2d<T, Coefficients>::hooke(const plane_t plane) const noexcept {
    const T E = this->E(plane);
    const T nu = this->nu(plane);
    return {     E / (T{1} - nu*nu), 
            nu * E / (T{1} - nu*nu),
        T{0.5} * E / (T{1} + nu) };
}
};

#endif