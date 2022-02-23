#ifndef NONLOCAL_SOLVERS_BOUNDARY_CONDITION_HPP
#define NONLOCAL_SOLVERS_BOUNDARY_CONDITION_HPP

#include <array>
#include <string>
#include <functional>

namespace nonlocal {

enum class boundary_type : uint8_t { FIRST_KIND, SECOND_KIND };

template<class T, class B, size_t DoF>
struct boundary_condition final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

    struct boundary_pair final {
        B type = B(boundary_type::SECOND_KIND);
        std::function<T(const std::array<T, 2>&)> func = [](const std::array<T, 2>&) constexpr noexcept { return 0; };
    };

    std::array<boundary_pair, DoF> data;

    static constexpr size_t degrees_of_freedom() noexcept { return DoF; }

    B type(const size_t b) const { return data[b].type; }
    const std::function<T(const std::array<T, 2>&)>& func(const size_t b) const { return data[b].func; }

    bool contains_condition_second_kind() const {
        return std::any_of(data.cbegin(), data.cend(), [](const boundary_pair& pair) { return pair.type == B(boundary_type::SECOND_KIND); });
    }
};

}

#endif