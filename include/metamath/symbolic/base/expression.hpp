#ifndef SYMBOLIC_EXPRESSION_HPP
#define SYMBOLIC_EXPRESSION_HPP

#define SYMBOLIC_NAMESPACE metamath::symbolic

namespace SYMBOLIC_NAMESPACE {

template<class E>
struct expression {
    constexpr E& operator()() noexcept {
        return static_cast<E&>(*this);
    }

    constexpr const E& operator()() const noexcept {
        return static_cast<const E&>(*this);
    }
};

}

#endif