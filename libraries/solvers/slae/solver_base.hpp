#pragma once

#include <metamath/linear/linear.hpp>
#include <parallel/OMP_utils.hpp>
#include <parallel/MPI_utils.hpp>

#include <optional>

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
class solver_base {
    const metamath::linear::sparse_matrix<T, I, J>& _matrix;
    parallel::MPI_ranges _process_rows;
    size_t _threads_count = parallel::threads_count();

public:
    using entity_t = metamath::types::container_type_t<T>;
    using floating_point_t = metamath::types::container_type_t<entity_t>;
    static_assert(std::is_floating_point_v<floating_point_t>, "Solver base class requires floating point type for calculations.");

    explicit solver_base(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _matrix{matrix}
        , _process_rows{parallel::rows_distribution(matrix.rows())} {}
    virtual ~solver_base() noexcept = default;

    const metamath::linear::sparse_matrix<T, I, J>& matrix() const noexcept {
        return _matrix;
    }

    const parallel::MPI_ranges& processes_ranges() const noexcept {
        return _process_rows;
    }

    std::ranges::iota_view<size_t, size_t> process_rows(const size_t process = parallel::MPI_rank()) const {
        return _process_rows.get(process);
    }

    size_t threads_count() const noexcept {
        return _threads_count;
    }

    virtual void set_threads_count(const size_t threads_count) {
        _threads_count = threads_count;
    }

    virtual std::vector<entity_t> solve(
        const std::vector<entity_t>& b,
        const std::optional<std::vector<entity_t>>& x0 = std::nullopt) const = 0;
};

}