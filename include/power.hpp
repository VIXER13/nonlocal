#ifndef POWER_HPP
#define POWER_HPP

namespace math_meta {

template<int64_t N, class Type> constexpr typename std::enable_if<(N > 0 && N % 2 == 0), Type>::type power(Type);
template<int64_t N, class Type> constexpr typename std::enable_if<(N > 0 && N % 2 == 1), Type>::type power(Type);

template<int64_t N, class Type>
constexpr typename std::enable_if<(N < 0), Type>::type power(Type x) {
    return 1 / power<-N>(x);
}

template<int64_t N, class Type>
constexpr typename std::enable_if<(N == 0), Type>::type power(Type) {
    return 1;
}

template<int64_t N, class Type>
constexpr typename std::enable_if<(N > 0 && N % 2 == 0), Type>::type power(Type x) {
    Type temp = power<N / 2>(x);
    return temp * temp;
}

template<int64_t N, class Type>
constexpr typename std::enable_if<(N > 0 && N % 2 == 1), Type>::type power(Type x) {
    return power<N - 1>(x) * x;
}

};

#endif