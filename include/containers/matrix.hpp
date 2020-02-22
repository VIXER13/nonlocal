#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

template<class Type, class Alloc = std::allocator<Type>>
class matrix {
    std::vector<Type, Alloc> matr;
    size_t row = 0, col = 0;

public:
    explicit matrix() noexcept = default;

    explicit matrix(const size_t rows, const size_t cols) :
        matr(rows * cols), row(cols ? rows : 0), col(rows ? cols : 0) {}

    explicit matrix(const size_t rows, const size_t cols, const Type value) :
        matr(rows * cols, value), row(cols ? rows : 0), col(rows ? cols : 0) {}

    size_t rows() const noexcept { return row; }
    size_t cols() const noexcept { return col; }
    size_t size() const noexcept { return matr.size(); }

    void resize(const size_t rows, const size_t cols) {
        row = cols ? rows : 0;
        col = rows ? cols : 0;
        matr.resize(row * col);
    }

    void resize(const size_t rows, const size_t cols, const Type value) {
        row = cols ? rows : 0;
        col = rows ? cols : 0;
        matr.resize(row * col, value);
        matr.shrink_to_fit();
    }

    const Type& operator()(const size_t i, const size_t j) const noexcept { return matr[i*col+j]; }
          Type& operator()(const size_t i, const size_t j)       noexcept { return matr[i*col+j]; }
        
    const Type* data() const noexcept { return matr.data(); }
          Type* data()       noexcept { return matr.data(); }
};

template<class Type, class Alloc = std::allocator<Type>>
std::ostream& operator<<(std::ostream &os, const matrix<Type, Alloc> &matr) {
    for(size_t i = 0; i < matr.rows(); ++i) {
        for(size_t j = 0; j < matr.cols(); ++j)
            os << matr(i, j) << " ";
        if(i < matr.rows())
            os << std::endl;
    }
    return os;
}

#endif