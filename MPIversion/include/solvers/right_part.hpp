#ifndef NONLOCAL_SOLVERS_RIGHT_PART_HPP
#define NONLOCAL_SOLVERS_RIGHT_PART_HPP

#include <array>
#include <functional>

namespace nonlocal {

template<class T, size_t DoF>
struct right_part final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

    struct function final {
        std::function<T(const std::array<T, 2>&)> func = [](const std::array<T, 2>&) constexpr noexcept { return 0; };
        T operator()(const std::array<T, 2>& x) const { return func(x); }
    };

    std::array<function, DoF> func;

    const function& operator[](const size_t i) const { return func[i]; }
};

}

#endif