#pragma once

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <limits>

namespace nonlocal::unit_tests {

template<std::floating_point T>
std::shared_ptr<mesh::mesh_1d<T>> make_small_mesh() {
    return std::make_shared<mesh::mesh_1d<T>>(
        metamath::finite_element::make_element_1d_integrated<T>(1, 2),
        std::vector<mesh::segment_data<T>>{ { .length = T{ 2 }, .search_radius = T{ 2 }, .elements = 2 } });
}

template<class Vector, std::floating_point T, size_t N>
void check_values(const Vector& values, const std::array<T, N>& expected, const T epsilon = T{ 1e-12 }) {
    using namespace boost::ut;
    expect(eq(size_t(values.size()), N)) << fatal;
    for (const size_t i : std::ranges::iota_view{ 0zu, N })
        expect(approx(values[i], expected[i], epsilon)) << " at node " << i;
}

template<std::floating_point T, class Solve, class Read_Field, class Read_Derived, class Expected_Field, class Expected_Derived>
void check_mesh_convergence(const T length, const size_t quadrature_order, const Solve& solve,
                            const Read_Field& read_field, const Read_Derived& read_derived,
                            const Expected_Field& expected_field, const Expected_Derived& expected_derived) {
    using namespace boost::ut;
    std::array<T, 4> previous_errors{};
    std::fill(previous_errors.begin(), previous_errors.end(), std::numeric_limits<T>::max());
    for (const size_t elements : { 10zu, 20zu, 40zu }) {
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(
            metamath::finite_element::make_element_1d_integrated<T>(1, quadrature_order),
            std::vector<mesh::segment_data<T>>{ { .length = length, .elements = elements } });
        auto solution = std::invoke(solve, mesh);
        const auto& field = std::invoke(read_field, solution);
        const auto& derived = std::invoke(read_derived, solution);
        const auto field_expected = mesh::utils::discrete<T>(*mesh, expected_field);
        const auto derived_expected = mesh::utils::discrete<T>(*mesh, expected_derived);
        // Preserve the four existing discrete error metrics and their refinement checks.
        const auto errors = std::array{
            L2_norm(field, field_expected), L2_norm(derived, derived_expected),
            max_error(field, field_expected), max_error(derived, derived_expected)
        };
        for (const size_t i : std::ranges::iota_view{ 0zu, errors.size() })
            expect(lt(errors[i], previous_errors[i]));
        previous_errors = errors;
    }
}

} // namespace nonlocal::unit_tests
