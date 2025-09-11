#include "init_2d_mesh.hpp"

namespace unit_tests {

template<>
std::shared_ptr<nonlocal::mesh::mesh_2d<double, int64_t>> init_2d_mesh(std::stringstream& stream, const nonlocal::mesh::mesh_format format) {
    return std::make_shared<nonlocal::mesh::mesh_2d<double, int64_t>>(stream, format);
}

}
