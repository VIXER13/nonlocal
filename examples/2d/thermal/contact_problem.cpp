#include "thermal_parameters_2d.hpp"

#include <unordered_map>
#include <iostream>

namespace {

using namespace nonlocal;

using T = double;

}

int main() {
    thermal::parameters_2d<T> params = {
        {"A", {.model = {.local_weight = 0.1}}},
        {"B", {.model = {.local_weight = 0.2}}}
    };

    std::cout << params["B"].model.local_weight << std::endl;

    return 0;
}