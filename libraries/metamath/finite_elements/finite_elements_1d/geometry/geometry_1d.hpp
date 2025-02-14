#pragma once

#include "geometry_1d_base.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Shape_Type, auto... Args>
class geometry_1d : public geometry_1d_base<T>,
                    public Shape_Type<T, Args...> {
protected:
    using shape_t = Shape_Type<T, Args...>;
    static_assert(shape_t::boundary.size() == 2, "Wrong number of boundaries.");

public:
    ~geometry_1d() override = default;

    T boundary(const side_1d bound) const override { return shape_t::boundary[size_t(bound)]; }
};

}