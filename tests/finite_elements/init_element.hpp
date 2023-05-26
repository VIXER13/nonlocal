#include "metamath.hpp"

namespace unit_tests {

// Element initialization is separated into a separate translation unit to avoid the SIOF problem

template<class T>
std::array<std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>, 6> init_elements();

}