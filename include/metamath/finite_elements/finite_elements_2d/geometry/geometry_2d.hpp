#ifndef FINITE_ELEMENT_GEOMETRY_2D_HPP
#define FINITE_ELEMENT_GEOMETRY_2D_HPP

#include "geometry_2d_base.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Shape_Type, auto... Args>
class geometry_2d : public geometry_2d_base<T>,
                    public Shape_Type<T, Args...> {
    using shape_t = Shape_Type<T, Args...>;
    static_assert(shape_t::boundary.size() == 4, "Wrong number of boundaries.");

public:
    ~geometry_2d() override = default;

    T boundary(const side_2d bound, const T x) const override { return shape_t::boundary[size_t(bound)](x); }
};

}

#endif