#include "init_2d_mesh.hpp"

namespace {

template<std::floating_point T, std::signed_integral I>
std::shared_ptr<nonlocal::mesh::mesh_2d<T, I>> init(const std::string& path_to_mesh) {
    return std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(path_to_mesh);
}

}

namespace unit_tests {

template<>
std::shared_ptr<nonlocal::mesh::mesh_2d<double, int64_t>> init_2d_mesh(const std::string& path_to_mesh) {
    return init<double, int64_t>(path_to_mesh);
}

}
