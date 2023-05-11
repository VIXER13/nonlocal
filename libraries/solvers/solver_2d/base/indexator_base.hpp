#ifndef NONLOCAL_INDEXATOR_BASE_HPP
#define NONLOCAL_INDEXATOR_BASE_HPP

#include "mesh_runner_types.hpp"

#include <array>
#include <vector>

namespace nonlocal {

template<size_t DoF>
class indexator_base {
    std::array<std::array<std::vector<bool>, DoF>, 2> _flags;
    const bool _is_symmetric;

protected:
    template<class Callback>
    static void check_flag(std::vector<bool>& flags, const size_t col, const Callback& callback);

    explicit indexator_base(const size_t size, const bool is_symmetric);

    bool is_symmetric() const noexcept;

    std::array<std::vector<bool>, DoF>& flags(const matrix_part part) noexcept;

public:
    virtual ~indexator_base() noexcept = default;

    void reset(const size_t node);
};

template<size_t DoF>
template<class Callback>
void indexator_base<DoF>::check_flag(std::vector<bool>& flags, const size_t col, const Callback& callback) {
    if (!flags[col]) {
        callback();
        flags[col] = true;
    }
}

template<size_t DoF>
indexator_base<DoF>::indexator_base(const size_t size, const bool is_symmetric)
    : _is_symmetric{is_symmetric} {
    for(auto& part : _flags)
        for(auto& dof : part)
            dof.resize(size);
}

template<size_t DoF>
bool indexator_base<DoF>::is_symmetric() const noexcept {
    return _is_symmetric;
}

template<size_t DoF>
std::array<std::vector<bool>, DoF>& indexator_base<DoF>::flags(const matrix_part part) noexcept {
    return _flags[size_t(part)];
}

template<size_t DoF>
void indexator_base<DoF>::reset(const size_t node) {
    for(const size_t i : std::ranges::iota_view{0u, DoF}) {
        std::fill(flags(matrix_part::BOUND)[i].begin(), flags(matrix_part::BOUND)[i].end(), false);
        std::fill(std::next(flags(matrix_part::INNER)[i].begin(), is_symmetric() ? DoF * node : 0), flags(matrix_part::INNER)[i].end(), false);
    }
}

}

#endif