#pragma once

namespace metamath::types {

template<class... Ts>
struct visitor : Ts... {
    using Ts::operator()...;
};

}