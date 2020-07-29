#ifndef MESH_2D_HPP
#define MESH_2D_HPP

#include <memory>
#include <string>
#include <fstream>
#include <exception>
#include "metamath/metamath.hpp"

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

enum class vtk_element_number : uintmax_t {
    LINEAR = 3,
    QUADRATIC = 21,
    TRIANGLE = 5,
    QUADRATIC_TRIANGLE = 22,
    BILINEAR = 9,
    QUADRATIC_SERENDIPITY = 23,
    QUADRATIC_LAGRANGE = 28
};

template<class Type, class Index = int32_t>
class mesh_2d {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    static_assert(std::is_integral_v<Index>, "The Index must be integral.");

    std::vector<std::array<Type, 2>> _nodes;                  // Глобальные кооридинаты узлов.
    std::vector<std::vector<Index>> _elements;                // Массив с глобальными номерами узлов для каждого элемента
                                                              // _elements[e][i] --- номер i-го узла элемента под номером e.
    std::vector<std::vector<std::vector<Index>>> _boundaries; // Массив с группами элементов на границе. Считаем, что на гранях элементов,
                                                              // которые располагаются вблизи границы, заданы одномерные элементы, поэтому
                                                              // _boundaries[b][e][i] - обращение к b-ой группе, i-ому узлу элемента под номером e.
    std::vector<std::vector<element_1d_t>> _elements_1d_type; // Массив с индексами одномерных элементов.
    std::vector<element_2d_t> _elements_2d_type;              // Массив с индексами двумерных элементов.
    std::vector<std::vector<Index>> _elements_neighbors,      // Ближайшие соседи элементов.
                                    _nodes_neighbors;         // Ближайшие соседи узлов.

    template<class T, template<class> class Quadrature_Type>
    using quadrature = metamath::finite_element::quadrature<T, Quadrature_Type>;
    template<class T>
    using gauss1 = metamath::finite_element::gauss1<T>;
    template<class T>
    using gauss2 = metamath::finite_element::gauss2<T>;
    template<class T>
    using gauss3 = metamath::finite_element::gauss3<T>;

    template<class T, template<class> class Element_Type>
    using element_1d_integrate = metamath::finite_element::element_1d_integrate<T, Element_Type>;
    template<class T>
    using linear = metamath::finite_element::linear<T>;
    template<class T>
    using quadratic = metamath::finite_element::quadratic<T>;

    template<class T, template<class> class Element_Type>
    using element_2d_integrate = metamath::finite_element::element_2d_integrate<T, Element_Type>;
    template<class T>
    using triangle = metamath::finite_element::triangle<T>;
    template<class T>
    using quadratic_triangle = metamath::finite_element::quadratic_triangle<T>;
    template<class T>
    using bilinear = metamath::finite_element::bilinear<T>;
    template<class T>
    using quadratic_serendipity = metamath::finite_element::quadratic_serendipity<T>;
    template<class T>
    using quadratic_lagrange = metamath::finite_element::quadratic_lagrange<T>;

    template<size_t... I>
    void read_element(std::ifstream& mesh_file, std::vector<Index>& element);
    void read_su2(const std::string& path);

    template<size_t... I>
    void write_element(std::ofstream& mesh_file, std::vector<Index>& element);

public:
    using fe_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<Type>>;
    using fe_2d_ptr = std::unique_ptr<metamath::finite_element::element_2d_integrate_base<Type>>;

private:
    // Наиболее употребительные одномерные и двумерные элементы
    std::array<fe_1d_ptr, 2> _elements_1d = {
        std::make_unique<element_1d_integrate<Type,    linear>>(quadrature<Type, gauss1>{}),
        std::make_unique<element_1d_integrate<Type, quadratic>>(quadrature<Type, gauss2>{}),
    };

    std::array<fe_2d_ptr, 5> _elements_2d = {
        std::make_unique<element_2d_integrate<Type,              triangle>>(quadrature<Type, gauss1>{}),
        std::make_unique<element_2d_integrate<Type,    quadratic_triangle>>(quadrature<Type, gauss2>{}),
        std::make_unique<element_2d_integrate<Type,              bilinear>>(quadrature<Type, gauss2>{}),
        std::make_unique<element_2d_integrate<Type, quadratic_serendipity>>(quadrature<Type, gauss3>{}),
        std::make_unique<element_2d_integrate<Type,    quadratic_lagrange>>(quadrature<Type, gauss3>{}),
    };

public:
    explicit mesh_2d(const std::string& path);
    void read_from_file(const std::string& path);

    void clear() noexcept;
    void shrink_to_fit() noexcept;

    size_t nodes_count() const noexcept;
    const std::array<Type, 2>& node(const size_t node) const noexcept;

    size_t elements_count() const noexcept;
    Index node_number(const size_t element, const size_t node) const noexcept;

    size_t boundary_groups_count() const noexcept;
    size_t elements_count(const size_t boundary) const noexcept;
    Index node_number(const size_t boundary, const size_t element, const size_t node) const noexcept;

    element_1d_t element_1d_type(const size_t boundary, const size_t element) const noexcept;
    const fe_1d_ptr& element_1d(const element_1d_t type) const noexcept;

    element_2d_t element_2d_type(const size_t element) const noexcept;
    const fe_2d_ptr& element_2d(const element_2d_t type) const noexcept;

    const std::vector<Index>& element_neighbors(const size_t element) const noexcept;
    const std::vector<Index>& node_neighbors(const size_t node) const noexcept;

    void find_elements_neighbors(const Type r); // Ищет соседние элементы, центры которых попали в зону влияния
    void find_nodes_neighbors(const Type r);    // Ищет соседние узлы, которые попадают в зону влияния

    void save_as_vtk(std::ofstream& fout) const;
};

template<class Type, class Index>
mesh_2d<Type, Index>::mesh_2d(const std::string& path) {
    read_from_file(path);
}

template<class Type, class Index>
void mesh_2d<Type, Index>::read_from_file(const std::string& path) {
    if (path.substr(path.size()-4) == ".su2") {
        read_su2(path);
    } else {
        throw std::domain_error{"Read format is not .su2."};
    }
}

template<class Type, class Index>
void mesh_2d<Type, Index>::clear() noexcept {
    _nodes.clear();
    _elements.clear();
    _boundaries.clear();
    _elements_1d_type.clear();
    _elements_2d_type.clear();
    _elements_neighbors.clear();
    _nodes_neighbors.clear();
}

template<class Type, class Index>
void mesh_2d<Type, Index>::shrink_to_fit() noexcept {
    _nodes.shrink_to_fit();
    _elements.shrink_to_fit();
    _boundaries.shrink_to_fit();
    _elements_1d_type.shrink_to_fit();
    _elements_2d_type.shrink_to_fit();
    _elements_neighbors.shrink_to_fit();
    _nodes_neighbors.shrink_to_fit();
}

template<class Type, class Index>
size_t mesh_2d<Type, Index>::nodes_count() const noexcept {
    return _nodes.size();
}

template<class Type, class Index>
const std::array<Type, 2>& mesh_2d<Type, Index>::node(const size_t node) const noexcept {
    return _nodes[node];
}

template<class Type, class Index>
size_t mesh_2d<Type, Index>::elements_count() const noexcept {
    return _elements.size();
}

template<class Type, class Index>
Index mesh_2d<Type, Index>::node_number(const size_t element, const size_t node) const noexcept {
    return _elements[element][node];
}

template<class Type, class Index>
size_t mesh_2d<Type, Index>::boundary_groups_count() const noexcept {
    return _boundaries.size();
}

template<class Type, class Index>
size_t mesh_2d<Type, Index>::elements_count(const size_t boundary) const noexcept {
    return _boundaries[boundary].size();
}

template<class Type, class Index>
Index mesh_2d<Type, Index>::node_number(const size_t boundary, const size_t element, const size_t node) const noexcept {
    return _boundaries[boundary][element][node];
}

template<class Type, class Index>
element_1d_t mesh_2d<Type, Index>::element_1d_type(const size_t boundary, const size_t element) const noexcept {
    return _elements_1d_type[boundary][element];
}

template<class Type, class Index>
const typename mesh_2d<Type, Index>::fe_1d_ptr& mesh_2d<Type, Index>::element_1d(const element_1d_t type) const noexcept {
    return _elements_1d[size_t(type)];
}

template<class Type, class Index>
element_2d_t mesh_2d<Type, Index>::element_2d_type(const size_t element) const noexcept {
    return _elements_2d_type[element];
}

template<class Type, class Index>
const typename mesh_2d<Type, Index>::fe_2d_ptr& mesh_2d<Type, Index>::element_2d(const element_2d_t type) const noexcept {
    return _elements_2d[size_t(type)];
}

template<class Type, class Index>
const std::vector<Index>& mesh_2d<Type, Index>::element_neighbors(const size_t element) const noexcept {
    return _elements_neighbors[element];
}

template<class Type, class Index>
const std::vector<Index>& mesh_2d<Type, Index>::node_neighbors(const size_t node) const noexcept {
    return _nodes_neighbors[node];
}

}

#endif