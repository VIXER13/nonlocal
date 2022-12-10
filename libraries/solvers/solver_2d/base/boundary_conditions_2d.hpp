#ifndef NONLOCAL_BOUNDARY_CONDITIONS_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITIONS_2D_HPP

#include <array>
#include <memory>
#include <string>
#include <unordered_map>

namespace nonlocal {

template<class T>
class boundary_condition_2d {
protected:
    constexpr explicit boundary_condition_2d() noexcept = default;

protected:
    static constexpr auto function_from_value(const T value) noexcept {
        return [value](const std::array<T, 2>&) constexpr noexcept { return value; };
    }

public:
    virtual ~boundary_condition_2d() noexcept = default;
    virtual T operator()(const std::array<T, 2>& x) const = 0;
};

template<class T>
class first_kind_2d : public virtual boundary_condition_2d<T> {
protected:
    constexpr explicit first_kind_2d() noexcept = default;

public:
    ~first_kind_2d() noexcept override = default;
};

template<class T>
class second_kind_2d : public virtual boundary_condition_2d<T> {
protected:
    constexpr explicit second_kind_2d() noexcept = default;

public:
    ~second_kind_2d() noexcept override = default;
};

template<class T, size_t DoF>
using boundary_conditions_2d = std::conditional_t<DoF == 1,
    std::unique_ptr<boundary_condition_2d<T>>,
    std::array<std::unique_ptr<boundary_condition_2d<T>>, DoF>
>;

template<class T, size_t DoF>
using boundaries_conditions_map_2d = std::unordered_map<std::string, boundary_conditions_2d<T, DoF>>;

}

#endif