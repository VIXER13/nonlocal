#ifndef MESH_UTILS_HPP
#define MESH_UTILS_HPP

// Данный модуль содержит функционал, который требует пересмотра

#include <cmath>
#include "mesh_2d.hpp"

namespace mesh::utils {

template<class T, size_t N>
T distance(const std::array<T, N>& A, const std::array<T, N>& B) noexcept {
    T sum = 0;
    for(size_t i = 0; i < N; ++i)
        sum += metamath::power<2>(A[i] - B[i]);
    return sqrt(sum);
}

bool is_trinagle(const element_2d_t type) noexcept {
    return type == element_2d_t::TRIANGLE || type == element_2d_t::QUADRATIC_TRIANGLE;
}

template<class T, size_t N>
std::array<T, N> operator*(std::array<T, N> arr, const T val) noexcept {
    for(size_t i = 0; i < N; ++i)
        arr[i] *= val;
    return val;
}

template<class T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& lhs, const std::array<T, N>& rhs) noexcept {
    for(size_t i = 0; i < N; ++i)
        lhs[i] += rhs[i];
    return lhs;
}

template<class T, class Index>
std::vector<std::array<T, 2>> get_centres_of_elements(const mesh_2d<T, Index>& mesh) {
    std::vector<std::array<T, 2>> centres(mesh.elements_count(), {});
    for(size_t element = 0; element < centres.size(); ++element) {
        const T x0 = is_trinagle(mesh.element_2d_type(element)) ? 1./3. : 0.;
        const auto& e = mesh.element_2d(mesh.element_2d_type(element));
        for(size_t node = 0; node < e->nodes_count(); ++node)
            centres[element] += mesh.node(node) * e->N(node, {x0, x0});
    }
    return std::move(centres);
}

}

#endif