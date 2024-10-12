#pragma once

#include <cstddef>
#include <vector>

namespace metamath::types {

template<class T>
class matrix : protected std::vector<T> {
    using _base = std::vector<T>;

    size_t _cols = 0;

public:
    using _base::size;
    using _base::begin;
    using _base::end;
    using _base::data;
    using _base::swap;
    using _base::operator[];

    explicit matrix(const size_t rows, const size_t cols);
    explicit matrix(const size_t rows, const size_t cols, const T& init);

    size_t rows() const noexcept;
    size_t cols() const noexcept;

    T& operator()(const size_t row, const size_t col);
    const T& operator()(const size_t row, const size_t col) const;
};

template<class T>
matrix<T>::matrix(const size_t rows, const size_t cols)
    : _base(rows * cols)
    , _cols{rows ? cols : 0} {}

template<class T>
matrix<T>::matrix(const size_t rows, const size_t cols, const T& init)
    : _base(rows * cols, init)
    , _cols{rows ? cols : 0} {}

template<class T>
size_t matrix<T>::rows() const noexcept {
    return size() / cols();
}

template<class T>
size_t matrix<T>::cols() const noexcept {
    return _cols;
}

template<class T>
T& matrix<T>::operator()(const size_t row, const size_t col) {
    return (*this)[row * cols() + col];
}

template<class T>
const T& matrix<T>::operator()(const size_t row, const size_t col) const {
    return (*this)[row * cols() + col];
}

}