#pragma once

#include "traits.hpp"

#include <concepts>
#include <memory>

namespace metamath::types {

template<class T>
requires copyable<T>
class copyable_uptr {
    std::unique_ptr<T> _ptr{};

public:
    constexpr copyable_uptr() noexcept = default;

    constexpr explicit copyable_uptr(std::unique_ptr<T>&& ptr) noexcept
        : _ptr(std::move(ptr)) {}

    constexpr copyable_uptr(const copyable_uptr& other) noexcept
        : _ptr(other._ptr ? other._ptr->copy() : nullptr) {}

    constexpr copyable_uptr(copyable_uptr&& other) noexcept
        : _ptr(std::move(other._ptr)) {}

    T* operator->() noexcept { return _ptr.get(); }
    const T* operator->() const noexcept { return _ptr.get(); }
    T& operator*() noexcept { return *_ptr; }
    const T& operator*() const noexcept { return *_ptr; }
    void reset(std::unique_ptr<T>&& ptr = nullptr) noexcept { _ptr.reset(ptr.release()); }
};

} // namespace metamath::types
