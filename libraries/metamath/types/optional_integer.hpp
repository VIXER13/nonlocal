#include <concepts>
#include <cstdint>
#include <climits>
#include <optional>
#include <stdexcept>

namespace metamath::types {

template<std::integral T>
class optional_integer final {
    bool _is_init : 1 = false;
    T    _value   : CHAR_BIT * sizeof(T) - 1 = T{0};

public:
    constexpr optional_integer() noexcept = default;
    constexpr optional_integer(const std::nullopt_t) noexcept {}

    template<class U>
    constexpr optional_integer(U&& other) noexcept
        : _is_init{true}
        , _value(other) {}

    constexpr optional_integer<T>& operator=(const std::nullopt_t) noexcept {
        _is_init = false;
        _value = 0;
        return *this;
    }

    template<class U>
    constexpr optional_integer<T>& operator=(U&& other) noexcept {
        _is_init = true;
        _value = other;
        return *this;
    }

    constexpr explicit operator bool() const noexcept {
        return _is_init;
    }

    constexpr operator T() const {
        if (!_is_init)
            throw std::runtime_error{"Attempting to access an uninitialized value."};
        return _value;
    }
};

}