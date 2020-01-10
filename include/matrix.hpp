// В данном модуле реализован класс матриц.
// Данный класс часто бывает удобен для некоторых мелких нужд, 
// ради которых не хочется подключать полноценную библиотеку.
// Он не претендует на полноту и оптимальность.

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

template<class Type, class Alloc = std::allocator<Type>>
class matrix
{
public:
    matrix() :
        row(0), col(0), matr(0) {}

    matrix(const size_t rows, const size_t cols) :
        row(cols == 0 ? 0 : rows), col(rows == 0 ? 0 : cols), matr(rows * cols) {}

    matrix(const size_t rows, const size_t cols, const Type value) :
        row(cols == 0 ? 0 : rows), col(rows == 0 ? 0 : cols), matr(rows * cols, value) {}

    matrix(const matrix<Type, Alloc> &other) :
        row(other.row), col(other.col), matr(other.matr) {}

    matrix(matrix<Type, Alloc> &&other) :
        row(other.row), col(other.col), matr(std::move(other.matr)) {}

    matrix<Type, Alloc>& operator=(const matrix<Type, Alloc> &other)
    {
        row = other.row;
        col = other.col;
        matr = other.matr;
        return *this;
    }

    matrix<Type, Alloc>& operator=(matrix<Type, Alloc> &&other)
    {
        row = other.row;
        col = other.col;
        matr = other.matr;
        return *this;
    }

    size_t rows() const { return row; }
    size_t cols() const { return col; }
    size_t size() const { return matr.size(); }

    void reserve(const size_t capacity)
    {
        if(capacity > row * col)
            matr.reserve(capacity);
    }

    void resize(const size_t rows, const size_t cols)
    {
        row = cols == 0 ? 0 : rows;
        col = rows == 0 ? 0 : cols;
        matr.resize(row * col);
    }

    void resize(const size_t rows, const size_t cols, const Type value)
    {
        row = cols == 0 ? 0 : rows;
        col = rows == 0 ? 0 : cols;
        matr.resize(row * col, value);
    }

    Type  operator()(const size_t i, const size_t j) const { return matr[i*col+j]; }
    Type& operator()(const size_t i, const size_t j)       { return matr[i*col+j]; }
        
    const Type* data() const { return matr.data(); }
          Type* data()       { return matr.data(); }

    void clear()
    {
        row = col = 0;
        matr.clear();
    }

private:
    size_t row, col;
    std::vector<Type, Alloc> matr;
};

template<class Type, class Alloc = std::allocator<Type>>
std::ostream& operator<<(std::ostream &os, const matrix<Type, Alloc> &matr)
{
    for(size_t i = 0; i < matr.rows(); ++i)
    {
        for(size_t j = 0; j < matr.cols(); ++j)
            os << matr(i, j) << " ";
        os << std::endl;
    }
    return os;
}

#endif