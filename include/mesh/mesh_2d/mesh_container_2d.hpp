#ifndef MESH_CONTAINER_2D_HPP
#define MESH_CONTAINER_2D_HPP

#include "metamath.hpp"
#include <memory>
#include <string>
#include <fstream>
#include <exception>
#include <sstream>
#include <iostream>

namespace nonlocal::mesh {

// Наиболее употребительные одномерные элементы
enum class element_1d_t : uint8_t {
    LINEAR,
    QUADRATIC
};

// Наиболее употребительные двумерные элементы
enum class element_2d_t : uint8_t {
    TRIANGLE,
    QUADRATIC_TRIANGLE,
    BILINEAR,
    QUADRATIC_SERENDIPITY,
    QUADRATIC_LAGRANGE
};

bool is_trinagle(const element_2d_t type) noexcept {
    return type == element_2d_t::TRIANGLE || type == element_2d_t::QUADRATIC_TRIANGLE;
}

enum class vtk_element_number : uintmax_t {
    LINEAR = 3,
    QUADRATIC = 21,
    TRIANGLE = 5,
    QUADRATIC_TRIANGLE = 22,
    BILINEAR = 9,
    QUADRATIC_SERENDIPITY = 23,
    QUADRATIC_LAGRANGE = 28
};

template<class T, class I = int32_t>
class mesh_2d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    static_assert(std::is_integral_v<I>, "The I must be integral.");

    std::vector<std::array<T, 2>> _nodes;
    std::vector<std::vector<I>>   _elements;
    std::vector<element_2d_t>     _elements_2d_type;

    std::vector<std::string>                                     _boundaries_names;
    std::unordered_map<std::string, std::vector<std::vector<I>>> _boundaries;
    std::unordered_map<std::string, std::vector<element_1d_t>>   _elements_1d_type;

    template<size_t... K, class Stream>
    void read_element(Stream& mesh_file, std::vector<I>& element);
    template<class Stream>
    void read_su2(Stream& mesh_file);

    template<size_t K0, size_t... K>
    void write_element(std::ofstream& mesh_file, const std::vector<I>& element) const;

public:
    using Finite_Element_1D_Ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;
    using Finite_Element_2D_Ptr = std::unique_ptr<metamath::finite_element::element_2d_integrate_base<T>>;

private:
    template<class U, template<class, auto...> class Quadrature_Type, auto... Args>
    using quadrature = metamath::finite_element::quadrature_1d<U, Quadrature_Type, Args...>;
    template<class U, size_t N>
    using gauss = metamath::finite_element::gauss<U, N>;

    template<class U, template<class, auto...> class Element_Type, auto... Args>
    using element_1d_integrate = metamath::finite_element::element_1d_integrate<U, Element_Type, Args...>;
    template<class U, size_t N>
    using lagrangian_element_1d = metamath::finite_element::lagrangian_element_1d<U, N>;

    template<class U, template<class, auto...> class Element_Type, auto... Args>
    using element_2d_integrate = metamath::finite_element::element_2d_integrate<U, Element_Type, Args...>;
    template<class U, size_t Order>
    using triangle = metamath::finite_element::triangle<U, Order>;
    template<class U, size_t Order>
    using serendipity = metamath::finite_element::serendipity<U, Order>;
    template<class U, size_t N, size_t M>
    using lagrangian_element_2d = metamath::finite_element::lagrangian_element_2d<U, N, M>;

    static std::array<Finite_Element_1D_Ptr, 2> make_default_1d_elements() {
        return { std::make_unique<element_1d_integrate<T, lagrangian_element_1d, 1>>(quadrature<T, gauss, 1>{}),
                 std::make_unique<element_1d_integrate<T, lagrangian_element_1d, 2>>(quadrature<T, gauss, 2>{}) };
    }

    static std::array<Finite_Element_2D_Ptr, 5> make_default_2d_elements() {
        return { std::make_unique<element_2d_integrate<T, triangle, 1>>(quadrature<T, gauss, 1>{}),
                 std::make_unique<element_2d_integrate<T, triangle, 2>>(quadrature<T, gauss, 2>{}),
                 std::make_unique<element_2d_integrate<T, serendipity, 1>>(quadrature<T, gauss, 2>{}),
                 std::make_unique<element_2d_integrate<T, serendipity, 2>>(quadrature<T, gauss, 3>{}),
                 std::make_unique<element_2d_integrate<T, lagrangian_element_2d, 2, 2>>(quadrature<T, gauss, 3>{}) };
    }

    // Наиболее употребительные одномерные и двумерные элементы
    std::array<Finite_Element_1D_Ptr, 2> _elements_1d = make_default_1d_elements();
    std::array<Finite_Element_2D_Ptr, 5> _elements_2d = make_default_2d_elements();

public:
    explicit mesh_2d(const std::string& path);
    mesh_2d(const mesh_2d& other);
    mesh_2d(mesh_2d&&) = default;

    void read_from_file(const std::string& path);

    size_t nodes_count() const noexcept;
    const std::array<T, 2>& node(const size_t node) const noexcept;
    const std::vector<std::array<T, 2>>& nodes() const noexcept;

    size_t elements_count() const noexcept;
    I node_number(const size_t element, const size_t node) const noexcept;
    size_t nodes_count(const size_t element) const noexcept;

    size_t boundary_groups_count() const noexcept;
    const std::vector<std::string>& boundary_names() const noexcept;
    size_t elements_count(const std::string& boundary) const;
    I node_number(const std::string& boundary, const size_t element, const size_t node) const noexcept;

    element_1d_t element_1d_type(const std::string& boundary, const size_t element) const noexcept;
    const Finite_Element_1D_Ptr& element_1d(const element_1d_t type) const noexcept;
    const Finite_Element_1D_Ptr& element_1d(const std::string& boundary, const size_t element) const noexcept;
    Finite_Element_1D_Ptr& element_1d(const element_1d_t type) noexcept;
    Finite_Element_1D_Ptr& element_1d(const std::string& boundary, const size_t element) noexcept;

    element_2d_t element_2d_type(const size_t element) const noexcept;
    const Finite_Element_2D_Ptr& element_2d(const element_2d_t type) const noexcept;
    const Finite_Element_2D_Ptr& element_2d(const size_t number) const noexcept;
    Finite_Element_2D_Ptr& element_2d(const element_2d_t type) noexcept;
    Finite_Element_2D_Ptr& element_2d(const size_t number) noexcept;

    void save_as_vtk(std::ofstream& mesh_file) const;
};

template<class T, class I>
mesh_2d<T, I>::mesh_2d(const std::string& path) {
    read_from_file(path);
}

template<class T, class I>
mesh_2d<T, I>::mesh_2d(const mesh_2d& other)
    : _nodes{other._nodes}
    , _elements{other._elements}
    , _elements_2d_type{other._elements_2d_type}
    , _boundaries_names{other._boundaries_names}
    , _boundaries{other._boundaries}
    , _elements_1d_type{other._elements_1d_type} {}

template<class T, class I>
void mesh_2d<T, I>::read_from_file(const std::string& path) {
    if (path.substr(path.size()-4) == ".su2") {
        std::ifstream file{path};
        read_su2(file);
    } else {
        throw std::domain_error{"Read format is not .su2."};
    }
}

template<class T, class I>
size_t mesh_2d<T, I>::nodes_count() const noexcept {
    return _nodes.size();
}

template<class T, class I>
const std::array<T, 2>& mesh_2d<T, I>::node(const size_t node) const noexcept {
    return _nodes[node];
}

template<class T, class I>
const std::vector<std::array<T, 2>>& mesh_2d<T, I>::nodes() const noexcept {
    return _nodes;
}

template<class T, class I>
size_t mesh_2d<T, I>::elements_count() const noexcept {
    return _elements.size();
}

template<class T, class I>
I mesh_2d<T, I>::node_number(const size_t element, const size_t node) const noexcept {
    return _elements[element][node];
}

template<class T, class I>
size_t mesh_2d<T, I>::nodes_count(const size_t element) const noexcept {
    return _elements[element].size();
}

template<class T, class I>
size_t mesh_2d<T, I>::boundary_groups_count() const noexcept {
    return _boundaries_names.size();
}

template<class T, class I>
const std::vector<std::string>& mesh_2d<T, I>::boundary_names() const noexcept {
    return _boundaries_names;
}

template<class T, class I>
size_t mesh_2d<T, I>::elements_count(const std::string& boundary) const {
    return _boundaries.at(boundary).size();
}

template<class T, class I>
I mesh_2d<T, I>::node_number(const std::string& boundary, const size_t element, const size_t node) const noexcept {
    return _boundaries.at(boundary)[element][node];
}

template<class T, class I>
element_1d_t mesh_2d<T, I>::element_1d_type(const std::string& boundary, const size_t element) const noexcept {
    return _elements_1d_type.at(boundary)[element];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const element_1d_t type) const noexcept {
    return _elements_1d[size_t(type)];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const std::string& boundary, const size_t element) const noexcept {
    return element_1d(element_1d_type(boundary, element));
}

template<class T, class I>
typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const element_1d_t type) noexcept {
    return _elements_1d[size_t(type)];
}

template<class T, class I>
typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const std::string& boundary, const size_t element) noexcept {
    return element_1d(element_1d_type(boundary, element));
}

template<class T, class I>
element_2d_t mesh_2d<T, I>::element_2d_type(const size_t element) const noexcept {
    return _elements_2d_type[element];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_2D_Ptr& mesh_2d<T, I>::element_2d(const element_2d_t type) const noexcept {
    return _elements_2d[size_t(type)];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_2D_Ptr& mesh_2d<T, I>::element_2d(const size_t number) const noexcept {
    return element_2d(element_2d_type(number));
}

template<class T, class I>
typename mesh_2d<T, I>::Finite_Element_2D_Ptr& mesh_2d<T, I>::element_2d(const element_2d_t type) noexcept {
    return _elements_2d[size_t(type)];
}

template<class T, class I>
typename mesh_2d<T, I>::Finite_Element_2D_Ptr& mesh_2d<T, I>::element_2d(const size_t number) noexcept {
    return element_2d(element_2d_type(number));
}

}

#endif