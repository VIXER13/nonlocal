#ifndef SYMDIFF_EXPRESSION_HPP
#define SYMDIFF_EXPRESSION_HPP

namespace metamath::symdiff {

template<class E>
struct expression {
    constexpr E& operator()() {
        return static_cast<E&>(*this);
    }

    constexpr const E& operator()() const {
        return static_cast<const E&>(*this);
    }
};

}

#endif