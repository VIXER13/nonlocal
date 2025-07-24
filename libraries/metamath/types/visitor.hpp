#pragma once

namespace metamath {

template<class... Ts>
struct visitor : Ts... {
    using Ts::operator()...;
};

}