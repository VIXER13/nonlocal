#pragma once

#include <array>
#include <functional>
#include <optional>

namespace nonlocal::solver_1d {

// A local scalar boundary law: value returns g(u), derivative returns its derivative dg/du.
// The stationary and transient assembly signs are defined in boundary_condition_nonlinear_1d.hpp.
template<class T>
struct nonlinear_boundary_operator_1d final {
    std::function<T(const T)> value;
    std::function<T(const T)> derivative;
};

// Robin (third-kind) condition for a local scalar operator:
//   a*du/dn + alpha*u = beta, where n is the outward normal and a is the material coefficient.
// Here alpha = robin_coefficient and beta = load.
// Heat convection: lambda*dT/dn + h*T = h*T_ambient.
//   alpha = h (heat transfer coefficient), beta = h*T_ambient (ambient temperature T_ambient).
// Mechanical spring: E*du/dn + k*u = p + k*u_ref.
//   alpha = k (spring stiffness), beta = p + k*u_ref (applied load p, reference displacement u_ref).
// These describe the linear boundary part; prescribed_value instead sets u directly.
template<class T>
struct boundary_data_1d final {
    std::optional<T> prescribed_value;
    T load = T{ 0 };
    std::optional<T> robin_coefficient;
    std::optional<nonlinear_boundary_operator_1d<T>> nonlinear;

    constexpr bool is_pure_load() const noexcept {
        return !prescribed_value && !robin_coefficient && !nonlinear;
    }
};

template<class T>
using boundaries_data_1d = std::array<boundary_data_1d<T>, 2>;

} // namespace nonlocal::solver_1d
