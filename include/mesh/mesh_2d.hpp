#ifndef MESH_2D_HPP
#define MESH_2D_HPP

#include <memory>
#include <string>
#include <fstream>
#include <exception>
#include <sstream>
#include "metamath.hpp"

namespace mesh {

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

    std::vector<std::array<T, 2>>            _nodes;              // Глобальные кооридинаты узлов.
    std::vector<std::vector<I>>              _elements;           // Массив с глобальными номерами узлов для каждого элемента
                                                                  // _elements[e][i] --- номер i-го узла элемента под номером e.
    std::vector<std::vector<std::vector<I>>> _boundaries;         // Массив с группами элементов на границе. Считаем, что на гранях элементов,
                                                                  // которые располагаются вблизи границы, заданы одномерные элементы, поэтому
                                                                  // _boundaries[b][e][i] - обращение к b-ой группе, i-ому узлу элемента под номером e.
    std::vector<std::vector<element_1d_t>>   _elements_1d_type;   // Массив с индексами одномерных элементов.
    std::vector<element_2d_t>                _elements_2d_type;   // Массив с индексами двумерных элементов.

    template<class U, template<class> class Quadrature_Type>
    using quadrature = metamath::finite_element::quadrature_1d<U, Quadrature_Type>;
    template<class U>
    using gauss1 = metamath::finite_element::gauss1<U>;
    template<class U>
    using gauss2 = metamath::finite_element::gauss2<U>;
    template<class U>
    using gauss3 = metamath::finite_element::gauss3<U>;

    template<class U, template<class> class Element_Type>
    using element_1d_integrate = metamath::finite_element::element_1d_integrate<U, Element_Type>;
    template<class U>
    using linear = metamath::finite_element::linear<U>;
    template<class U>
    using quadratic = metamath::finite_element::quadratic<U>;

    template<class U, template<class> class Element_Type>
    using element_2d_integrate = metamath::finite_element::element_2d_integrate<U, Element_Type>;
    template<class U>
    using triangle = metamath::finite_element::triangle<U>;
    template<class U>
    using quadratic_triangle = metamath::finite_element::quadratic_triangle<U>;
    template<class U>
    using bilinear = metamath::finite_element::bilinear<U>;
    template<class U>
    using quadratic_serendipity = metamath::finite_element::quadratic_serendipity<U>;
    template<class U>
    using quadratic_lagrange = metamath::finite_element::quadratic_lagrange<U>;

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
    static std::array<Finite_Element_1D_Ptr, 2> make_default_1d_elements() {
        return { std::make_unique<element_1d_integrate<T,    linear>>(quadrature<T, gauss1>{}),
                 std::make_unique<element_1d_integrate<T, quadratic>>(quadrature<T, gauss2>{}) };
    }

    static std::array<Finite_Element_2D_Ptr, 5> make_default_2d_elements() {
        return { std::make_unique<element_2d_integrate<T,              triangle>>(quadrature<T, gauss1>{}),
                 std::make_unique<element_2d_integrate<T,    quadratic_triangle>>(quadrature<T, gauss2>{}),
                 std::make_unique<element_2d_integrate<T,              bilinear>>(quadrature<T, gauss2>{}),
                 std::make_unique<element_2d_integrate<T, quadratic_serendipity>>(quadrature<T, gauss3>{}),
                 std::make_unique<element_2d_integrate<T,    quadratic_lagrange>>(quadrature<T, gauss3>{}) };
    }

    // Наиболее употребительные одномерные и двумерные элементы
    std::array<Finite_Element_1D_Ptr, 2> _elements_1d = make_default_1d_elements();
    std::array<Finite_Element_2D_Ptr, 5> _elements_2d = make_default_2d_elements();

    bool check_intersect(const std::array<T, 2>& a, const std::array<T, 2>& b) const;

public:
    explicit mesh_2d(const std::string& path);
    mesh_2d(const mesh_2d& other);
    mesh_2d(mesh_2d&&) = default;

    void read_from_file(const std::string& path);

    void clear() noexcept;
    void shrink_to_fit() noexcept;

    size_t nodes_count() const noexcept;
    const std::array<T, 2>& node(const size_t node) const noexcept;

    size_t elements_count() const noexcept;
    I node_number(const size_t element, const size_t node) const noexcept;
    size_t nodes_count(const size_t element) const noexcept;

    size_t boundary_groups_count() const noexcept;
    size_t elements_count(const size_t boundary) const noexcept;
    I node_number(const size_t boundary, const size_t element, const size_t node) const noexcept;

    element_1d_t element_1d_type(const size_t boundary, const size_t element) const noexcept;
    const Finite_Element_1D_Ptr& element_1d(const element_1d_t type) const noexcept;
    const Finite_Element_1D_Ptr& element_1d(const size_t boundary, const size_t element) const noexcept;

    element_2d_t element_2d_type(const size_t element) const noexcept;
    const Finite_Element_2D_Ptr& element_2d(const element_2d_t type) const noexcept;
    const Finite_Element_2D_Ptr& element_2d(const size_t number) const noexcept;

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
    , _boundaries{other._boundaries}
    , _elements_1d_type{other._elements_1d_type}
    , _elements_2d_type{other._elements_2d_type} {}

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
void mesh_2d<T, I>::clear() noexcept {
    _nodes.clear();
    _elements.clear();
    _boundaries.clear();
    _elements_1d_type.clear();
    _elements_2d_type.clear();
}

template<class T, class I>
void mesh_2d<T, I>::shrink_to_fit() noexcept {
    _nodes.shrink_to_fit();
    _elements.shrink_to_fit();
    _boundaries.shrink_to_fit();
    _elements_1d_type.shrink_to_fit();
    _elements_2d_type.shrink_to_fit();
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
    return _boundaries.size();
}

template<class T, class I>
size_t mesh_2d<T, I>::elements_count(const size_t boundary) const noexcept {
    return _boundaries[boundary].size();
}

template<class T, class I>
I mesh_2d<T, I>::node_number(const size_t boundary, const size_t element, const size_t node) const noexcept {
    return _boundaries[boundary][element][node];
}

template<class T, class I>
element_1d_t mesh_2d<T, I>::element_1d_type(const size_t boundary, const size_t element) const noexcept {
    return _elements_1d_type[boundary][element];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const element_1d_t type) const noexcept {
    return _elements_1d[size_t(type)];
}

template<class T, class I>
const typename mesh_2d<T, I>::Finite_Element_1D_Ptr& mesh_2d<T, I>::element_1d(const size_t boundary, const size_t element) const noexcept {
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

}

#endif