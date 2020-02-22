#ifndef ROWS_DIFFERENT_SIZES_HPP
#define ROWS_DIFFERENT_SIZES_HPP

#include <array>
#include <vector>

template<class Type, class Index = uint32_t, class Alloc_Data = std::allocator<Type>, class Alloc_Index = std::allocator<Index>>
class rows_different_sizes
{
    static_assert(std::is_integral<Index>::value, "The Index must be integral type.");

    std::vector<Type, Alloc_Data> matr;
    std::vector<Index, Alloc_Index> shift{0};

public:
    explicit rows_different_sizes() = default;

    size_t size() const noexcept { return matr.size(); }
    size_t rows() const noexcept { return shift.size()-1; }
    Index  cols(const size_t row) const noexcept { return shift[row+1] - shift[row]; }

    const Type& operator()(const size_t row, const size_t col) const noexcept { return matr[shift[row]+col]; }
          Type& operator()(const size_t row, const size_t col)       noexcept { return matr[shift[row]+col]; }

    const Type& back() const noexcept { return matr.back(); }
          Type& back()       noexcept { return matr.back(); }

    auto cbegin() const noexcept { return matr.cbegin(); }
    auto cend()   const noexcept { return matr.cend();   }

    void reserve(const size_t size, const size_t factor = 1) {
        matr.reserve(size * (factor ? factor : 1));
        shift.reserve(size);
    }

    template<size_t Size>
    void push_back(const std::array<Type, Size> &push) {
        if(matr.capacity() < matr.size() + Size)
            matr.reserve(matr.size() + Size);
        for(size_t i = 0; i < Size; ++i)
            matr.push_back(push[i]);
        shift.push_back(shift.back() + Size);
    }

    void push_back(const std::vector<Type> &push) {
        if(matr.capacity() < matr.size() + push.size())
            matr.reserve(matr.size() + push.size());
        for(size_t i = 0; i < push.size(); ++i)
            matr.push_back(push[i]);
        shift.push_back(shift.back() + push.size());
    }

    void clear() noexcept {
        matr.clear();
        shift.clear();
        shift.resize(1, 0);
    }

    void shrink_to_fit() {
        matr.shrink_to_fit();
        shift.shrink_to_fit();
    }    
};

template<class Type, class Index>
std::ostream& operator<<(std::ostream &os, const rows_different_sizes<Type, Index> &matr) {
    for(size_t i = 0; i < matr.rows(); ++i) {
        for(Index j = 0; j < matr.cols(i); ++j)
            os << matr(i, j) << " ";
        if(i != matr.rows()-1) 
            os << std::endl;
    }
    return os;
}

#endif