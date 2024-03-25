#ifndef UNIT_TESTS_BOUNDARY_CONDITIONS_HPP
#define UNIT_TESTS_BOUNDARY_CONDITIONS_HPP

#include <boost/ut.hpp>
#include <tuple>
#include <filesystem>

#include "logger.hpp"
#include "metamath.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

using T = double;
using I = int64_t;
static constexpr T epsilon = T{1e-8};

namespace nonstat_1d_tests {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

template <std::floating_point T>
struct time_data final {
    T time_step = T{0};       
    T initial_time = T{0};
    uint64_t steps_count = 0; 
};

template<class T, size_t Order>
using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, Order>;
template<class T, size_t Order>
using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, Order>;

}

#endif
