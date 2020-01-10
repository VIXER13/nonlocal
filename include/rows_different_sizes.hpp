#ifndef ROWS_DIFFERENT_SIZES_HPP
#define ROWS_DIFFERENT_SIZES_HPP

// Данный класс предназначен для хранения строк разных длин в виде разреженной матрицы вида: **********
//                                                                                           ******
//                                                                                           ********
//                                                                                           ************
// Интерфейс во многом повторяет интерфейс std::vector.
// Для возможной экономии памяти, есть возможность выбрать тип для сдвигов, который естественно должен быть целочисленным.

#include <array>
#include <vector>

template<class Type, class Index = uint32_t>
class rows_different_sizes
{
    static_assert(std::is_integral<Index>::value, "Index must be integral.");

    std::vector<Type> matr;
    std::vector<Index> shift;

public:
    rows_different_sizes() :
        shift(1, 0) {}

    rows_different_sizes(const rows_different_sizes<Type, Index> &other) :
        matr(other.matr), shift(other.shift) {}

    rows_different_sizes(rows_different_sizes<Type, Index> &&other) :
        matr(std::move(other.matr)), shift(std::move(other.shift)) {}

    const Type& operator()(const size_t row, const size_t col) const { return matr[shift[row]+col]; }
          Type& operator()(const size_t row, const size_t col)       { return matr[shift[row]+col]; }

    const Type& back() const { return matr.back(); }
          Type& back()       { return matr.back(); }

    auto cbegin() const { return matr.cbegin(); }
    auto cend()   const { return matr.cend();   }

    size_t size() const { return matr.size(); }
    size_t rows() const { return shift.size()-1; }
    Index  cols(const size_t row) const { return shift[row+1] - shift[row]; }

    void reserve(const size_t size, const size_t factor = 1)
    {
        matr.reserve(size * (factor ? factor : 1));
        shift.reserve(size);
    }

    template<size_t Size>
    void push_back(const std::array<Type, Size> &push)
    {
        if(matr.capacity() < matr.size() + Size)
            matr.reserve(matr.size() + Size);
        for(size_t i = 0; i < Size; ++i)
            matr.push_back(push[i]);
        shift.push_back(shift.back() + Size);
    }

    void push_back(const std::vector<Type> &push)
    {
        if(matr.capacity() < matr.size() + push.size())
            matr.reserve(matr.size() + push.size());
        for(size_t i = 0; i < push.size(); ++i)
            matr.push_back(push[i]);
        shift.push_back(shift.back() + push.size());
    }

    void clear()
    {
        matr.clear();
        shift.clear();
        shift.resize(1, 0);
    }

    void shrink_to_fit()
    {
        matr.shrink_to_fit();
        shift.shrink_to_fit();
    }    
};

template<class Type, class Index>
std::ostream& operator<<(std::ostream &os, const rows_different_sizes<Type, Index> &matr)
{
    for(size_t i = 0; i < matr.rows(); ++i)
    {
        for(Index j = 0; j < matr.cols(i); ++j)
            os << matr(i, j) << " ";
        if(i != matr.rows()-1) 
            os << std::endl;
    }
    return os;
}

#endif