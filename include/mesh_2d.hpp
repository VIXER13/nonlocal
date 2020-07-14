#ifndef MESH_2D_HPP
#define MESH_2D_HPP

// В данном модуле реализован простейший класс сеток с возможностью генерации прямоугольной сетки
// из билинейных, квадратичных серендиповых и кубических серендиповых элементов.
// Генерация сетки в дальнейшем будет убрана, сейчас же она встроена для тестирования и отладки программы.

#include <fstream>
#include "matrix.hpp"
#include "rows_different_sizes.hpp"
#include "metamath/power.hpp"
#include "metamath/finite_element/finite_element.hpp"

template<class Type, class Index = int32_t>
class mesh_2d
{
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");
    static_assert(std::is_integral_v<Index>, "The Index must be integral.");

    rows_different_sizes<Index> elements;               // Массив с глобальными номерами узлов для каждого элемента. elements(n, i) - обращение к i-му узлу, n-го элемента
    matrix<Type> nodes;                                 // Глобальные координаты узлов сетки. nodes(n, i) - обращение к i-ой компоненте n-го узла

    std::vector<rows_different_sizes<Index>> bounds;    // Группы границ, каждая из которых содержит массив номеров узлов для одномерных элементов
                                                        // bound[b](n, j) - обращение к b-ой группе, j-ому узлу n-го элемента в группе

    std::vector<uint8_t>              elements_2d_type; // Массив с индексами элементов.
    std::vector<std::vector<uint8_t>> elements_1d_type; // Массивы с индексами элементов на границах.

    std::vector<std::vector<Index>> neighbors;          // Номера соседних элементов. neighbors[el][i] - i-ый сосед el-го элемента.

    std::array<metamath::finite_element::element_1d_integrate_base<Type>*, 3> finite_elements_1d =
        { new metamath::finite_element::element_1d_integrate<Type, metamath::finite_element::linear>   ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss1>())),
          new metamath::finite_element::element_1d_integrate<Type, metamath::finite_element::quadratic>((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>())),
          new metamath::finite_element::element_1d_integrate<Type, metamath::finite_element::qubic>    ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>())) };

    std::array<metamath::finite_element::element_2d_integrate_base<Type>*, 6> finite_elements_2d =
        { new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::triangle>          ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss1>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss1>())),
          new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::quadratic_triangle>((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>())),
          new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::qubic_triangle>    ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss3>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss3>())),
          new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::bilinear>          ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss2>())),
          new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::quadratic_serendipity>((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss3>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss3>())),
          new metamath::finite_element::element_2d_integrate<Type, metamath::finite_element::qubic_serendipity>    ((metamath::finite_element::quadrature<Type, metamath::finite_element::gauss4>()), (metamath::finite_element::quadrature<Type, metamath::finite_element::gauss4>())) };

    void make_bilinear          (const Index Ex, const Index Ey, const Type Lx, const Type Ly);
    void make_quadratic_serendip(const Index Ex, const Index Ey, const Type Lx, const Type Ly);
    void make_qubic_serendip    (const Index Ex, const Index Ey, const Type Lx, const Type Ly);

public:
    enum elemType {BILINEAR, QUADRATIC_SERENDIP, QUBIC_SERENDIP};

    mesh_2d(const elemType type, const Index Ex, const Index Ey, const Type Lx, const Type Ly);
    mesh_2d(const std::string& path);

    void read_su2(const std::string& path);

    void clear();

    bool check_intersect_with_bound(const std::array<Type, 2>& A, const std::array<Type, 2>& B);
    void find_neighbors_for_elements(const Type r);
    void find_neighbors_for_nodes(const Type r);

    size_t nodes_count() const noexcept;
    size_t elements_count() const noexcept;
    size_t boundary_groups_count() const noexcept;

    Index node_number(const size_t el, const size_t i) const noexcept;                   // Глобальный номер узла под номером i элемента под номером el
    Type  coord(const size_t node, const size_t component) const noexcept;               // Значение компоненты узла под номером node

    const rows_different_sizes<Index>& boundary(const size_t b) const noexcept;          // Возвращает массив с номерами узлов на границе b

    uint8_t element_type(const size_t el) const noexcept;                                // возвращает номер типа элемента el
    const std::vector<uint8_t>& elements_on_bound_types(const size_t b) const noexcept;  // Возвращает массив с типами элементов на границе b

    const std::vector<Index>& neighbor(const size_t el) const noexcept;                  // Возвращает массив с номерами соседних элементов

          metamath::finite_element::element_1d_integrate_base<Type>* element_1d(const size_t el)       noexcept;
    const metamath::finite_element::element_1d_integrate_base<Type>* element_1d(const size_t el) const noexcept;
          metamath::finite_element::element_2d_integrate_base<Type>* element_2d(const size_t el)       noexcept;
    const metamath::finite_element::element_2d_integrate_base<Type>* element_2d(const size_t el) const noexcept;

    void print_to_file(const std::string& path, const Type *x) const;

    ~mesh_2d();    
};

template<class Type, class Index>
void mesh_2d<Type, Index>::make_bilinear(const Index Ex, const Index Ey, const Type Lx, const Type Ly)
{
    elements.reserve(Ex*Ey, 4);
    for(Index i = 0; i < Ey; ++i)
        for(Index j = 0; j < Ex; ++j)
            elements.push_back(std::array<Index, 4>({    i*(Ex+1)+j,
                                                         i*(Ex+1)+j+1,
                                                     (i+1)*(Ex+1)+j+1,
                                                     (i+1)*(Ex+1)+j   }));

    const Type hx = Lx/Ex, hy = Ly/Ey;
    nodes.resize((Ex+1)*(Ey+1), 2);
    for(Index i = 0; i < Ey+1; ++i)
        for(Index j = 0; j < Ex+1; ++j)
        {
            nodes(i*(Ex+1)+j, 0) = hx*j;
            nodes(i*(Ex+1)+j, 1) = hy*i;
        }

    bounds.resize(4);
    
    bounds[0].reserve(Ex, 2);
    for(Index i = 0; i < Ex; ++i)
        bounds[0].push_back(std::array<Index, 2>({ elements(i, 0), elements(i, 1) }));

    bounds[1].reserve(Ey, 2);
    for(Index i = Ex-1; i < Ex*Ey; i += Ex)
        bounds[1].push_back(std::array<Index, 2>({ elements(i, 1), elements(i, 2) }));

    bounds[2].reserve(Ex, 2);
    for(Index i = Ex*(Ey-1); i < Ex*Ey; ++i)
        bounds[2].push_back(std::array<Index, 2>({ elements(i, 3), elements(i, 2) }));

    bounds[3].reserve(Ey, 2);
    for(Index i = 0; i < Ex*(Ey-1)+1; i += Ex)
        bounds[3].push_back(std::array<Index, 2>({ elements(i, 0), elements(i, 3) }));

    elements_2d_type.resize(Ex*Ey);
    for(size_t i = 0; i < elements_2d_type.size(); ++i)
        elements_2d_type[i] = 3;

    elements_1d_type.resize(4);
    elements_1d_type[0].resize(Ex);
    elements_1d_type[1].resize(Ey);
    elements_1d_type[2].resize(Ex);
    elements_1d_type[3].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 0;

    //////////////////////////////////////////////////////////////////////////////////
    // КОСТЫЛЬ ДЛЯ ЗАДАЧИ НЕЙМАНА
    bounds.resize(6);

    bounds[4].reserve(Ex, 2);
    for(Index i = Ex * Ey / 2; i < Ex * Ey / 2 + Ex; ++i)
    bounds[4].push_back(std::array<Index, 2>{ elements(i, 0), elements(i, 1) });

    bounds[5].reserve(Ey, 2);
    for(Index i = Ex / 2 - 1; i < Ex * Ey; i += Ex)
        bounds[5].push_back(std::array<Index, 2>{ elements(i, 1), elements(i, 2) });


    elements_1d_type.resize(6);

    elements_1d_type[4].resize(Ex);
    elements_1d_type[5].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 0;
    //////////////////////////////////////////////////////////////////////////////////
}

template<class Type, class Index>
void mesh_2d<Type, Index>::make_quadratic_serendip(const Index Ex, const Index Ey, const Type Lx, const Type Ly)
{
    elements.reserve(Ex*Ey, 8);
    for(Index i = 0; i < Ey; ++i)
        for(Index j = 0; j < Ex; ++j)
            elements.push_back(std::array<Index, 8>({ (3*Ex+2)*i     + 2*j,
                                                      (3*Ex+2)*i     + 2*j+1,
                                                      (3*Ex+2)*i     + 2*j+2,
                                                      (3*Ex+2)*i     + 2*Ex+1 + j+1,
                                                      (3*Ex+2)*(i+1) + 2*j+2,
                                                      (3*Ex+2)*(i+1) + 2*j+1,
                                                      (3*Ex+2)*(i+1) + 2*j,
                                                      (3*Ex+2)*i     + 2*Ex+1 + j }));

    const Type hx = Lx/Ex, hy = Ly/Ey;
    nodes.resize((3*Ex+2)*Ey+2*Ex+1, 2);
    for(Index i = 0; i < Ey+1; ++i)
    {
        for(Index j = 0; j < 2*Ex+1; ++j)
        {
            nodes(i*(3*Ex+2)+j, 0) = 0.5*hx*j;
            nodes(i*(3*Ex+2)+j, 1) = hy*i;
        }

        if(i < Ey)
            for(Index j = 0; j < Ex+1; ++j)
            {
                nodes(i*(3*Ex+2)+j+2*Ex+1, 0) = hx*j;
                nodes(i*(3*Ex+2)+j+2*Ex+1, 1) = hy*i+0.5*hy;
            }
    }

    bounds.resize(4);
    
    bounds[0].reserve(Ex, 3);
    for(Index i = 0; i < Ex; ++i)
        bounds[0].push_back(std::array<Index, 3>({ elements(i, 0), elements(i, 1), elements(i, 2) }));

    bounds[1].reserve(Ey, 3);
    for(Index i = Ex-1; i < Ex*Ey; i += Ex)
        bounds[1].push_back(std::array<Index, 3>({ elements(i, 2), elements(i, 3), elements(i, 4) }));

    bounds[2].reserve(Ex, 3);
    for(Index i = Ex*(Ey-1); i < Ex*Ey; ++i)
        bounds[2].push_back(std::array<Index, 3>({ elements(i, 6), elements(i, 5), elements(i, 4) }));

    bounds[3].reserve(Ey, 3);
    for(Index i = 0; i < Ex*(Ey-1)+1; i += Ex)
        bounds[3].push_back(std::array<Index, 3>({ elements(i, 0), elements(i, 7), elements(i, 6) }));

    elements_2d_type.resize(Ex*Ey);
    for(size_t i = 0; i < elements_2d_type.size(); ++i)
        elements_2d_type[i] = 4;

    elements_1d_type.resize(4);
    elements_1d_type[0].resize(Ex);
    elements_1d_type[1].resize(Ey);
    elements_1d_type[2].resize(Ex);
    elements_1d_type[3].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 1;

    //////////////////////////////////////////////////////////////////////////////////
    // КОСТЫЛЬ ДЛЯ ЗАДАЧИ НЕЙМАНА
    bounds.resize(6);

    bounds[4].reserve(Ex, 2);
    for(Index i = Ex * Ey / 2; i < Ex * Ey / 2 + Ex; ++i)
    bounds[4].push_back(std::array<Index, 3>{ elements(i, 0), elements(i, 1), elements(i, 2) });

    bounds[5].reserve(Ey, 2);
    for(Index i = Ex / 2 - 1; i < Ex * Ey; i += Ex)
        bounds[5].push_back(std::array<Index, 3>{ elements(i, 2), elements(i, 3), elements(i, 4) });


    elements_1d_type.resize(6);

    elements_1d_type[4].resize(Ex);
    elements_1d_type[5].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 1;
    //////////////////////////////////////////////////////////////////////////////////
}

template<class Type, class Index>
void mesh_2d<Type, Index>::make_qubic_serendip(const Index Ex, const Index Ey, const Type Lx, const Type Ly)
{
    elements.reserve(Ex*Ey, 12);
    for(Index i = 0; i < Ey; ++i)
        for(Index j = 0; j < Ex; ++j)
            elements.push_back(std::array<Index, 12>({ (5*Ex+3)*i + 3*j,
                                                       (5*Ex+3)*i + 3*j+1,
                                                       (5*Ex+3)*i + 3*j+2,
                                                       (5*Ex+3)*i + 3*j+3,
                                                       (5*Ex+3)*i + 3*Ex+1 + j+1,
                                                       (5*Ex+3)*i + 4*Ex+2 + j+1,
                                                       (5*Ex+3)*(i+1) + 3*j+3,
                                                       (5*Ex+3)*(i+1) + 3*j+2,
                                                       (5*Ex+3)*(i+1) + 3*j+1,
                                                       (5*Ex+3)*(i+1) + 3*j,
                                                       (5*Ex+3)*i + 4*Ex+2 + j,
                                                       (5*Ex+3)*i + 3*Ex+1 + j }));

    const Type hx = Lx/Ex, hy = Ly/Ey;
    nodes.resize((5*Ex+3)*Ey+3*Ex+1, 2);
    for(Index i = 0; i < Ey+1; ++i)
    {
        for(Index j = 0; j < 3*Ex+1; ++j)
        {
			nodes(i*(5*Ex+3)+j, 0) = 1.0/3.0*hx*j;
			nodes(i*(5*Ex+3)+j, 1) = hy*i;
		}

        if(i < Ey)
        {
			for(Index j = 0; j < Ex+1; ++j)
            {
				nodes(i*(5*Ex+3)+j+3*Ex+1, 0) = hx*j;
				nodes(i*(5*Ex+3)+j+3*Ex+1, 1) = hy*i+1.0/3.0*hy;
			}

			for(Index j = 0; j < Ex+1; ++j)
            {
				nodes(i*(5*Ex+3)+j+4*Ex+2, 0) = hx*j;
				nodes(i*(5*Ex+3)+j+4*Ex+2, 1) = hy*i+2.0/3.0*hy;
			}
		}
    }

    bounds.resize(4);
    
    bounds[0].reserve(Ex, 4);
    for(Index i = 0; i < Ex; ++i)
        bounds[0].push_back(std::array<Index, 4>({ elements(i, 0), elements(i, 1), elements(i, 2), elements(i, 3) }));

    bounds[1].reserve(Ey, 4);
    for(Index i = Ex-1; i < Ex*Ey; i += Ex)
        bounds[1].push_back(std::array<Index, 4>({ elements(i, 3), elements(i, 4), elements(i, 5), elements(i, 6) }));

    bounds[2].reserve(Ex, 4);
    for(Index i = Ex*(Ey-1); i < Ex*Ey; ++i)
        bounds[2].push_back(std::array<Index, 4>({ elements(i, 9), elements(i, 8), elements(i, 7), elements(i, 6) }));

    bounds[3].reserve(Ey, 4);
    for(Index i = 0; i < Ex*(Ey-1)+1; i += Ex)
        bounds[3].push_back(std::array<Index, 4>({ elements(i, 0), elements(i, 11), elements(i, 10), elements(i, 9) }));

    elements_2d_type.resize(Ex*Ey);
    for(size_t i = 0; i < elements_2d_type.size(); ++i)
        elements_2d_type[i] = 5;

    elements_1d_type.resize(4);
    elements_1d_type[0].resize(Ex);
    elements_1d_type[1].resize(Ey);
    elements_1d_type[2].resize(Ex);
    elements_1d_type[3].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 2;

    //////////////////////////////////////////////////////////////////////////////////
    // КОСТЫЛЬ ДЛЯ ЗАДАЧИ НЕЙМАНА
    bounds.resize(6);

    bounds[4].reserve(Ex, 2);
    for(Index i = Ex * Ey / 2; i < Ex * Ey / 2 + Ex; ++i)
    bounds[4].push_back(std::array<Index, 4>{ elements(i, 0), elements(i, 1), elements(i, 2), elements(i, 3) });

    bounds[5].reserve(Ey, 2);
    for(Index i = Ex / 2 - 1; i < Ex * Ey; i += Ex)
        bounds[5].push_back(std::array<Index, 4>{ elements(i, 3), elements(i, 4), elements(i, 5), elements(i, 6) });


    elements_1d_type.resize(6);

    elements_1d_type[4].resize(Ex);
    elements_1d_type[5].resize(Ey);
    for(size_t i = 0; i < elements_1d_type.size(); ++i)
        for(size_t j = 0; j < elements_1d_type[i].size(); ++j)
            elements_1d_type[i][j] = 2;
    //////////////////////////////////////////////////////////////////////////////////
}

template<class Type, class Index>
mesh_2d<Type, Index>::mesh_2d(const elemType type, const Index Ex, const Index Ey, const Type Lx, const Type Ly) :
    neighbors(Ex*Ey)
{
    switch (type)
    {
    case elemType::BILINEAR:
        make_bilinear(Ex, Ey, Lx, Ly);
        break;

    case elemType::QUADRATIC_SERENDIP:
        make_quadratic_serendip(Ex, Ey, Lx, Ly);
        break;

    case elemType::QUBIC_SERENDIP:
        make_qubic_serendip(Ex, Ey, Lx, Ly);
        break;
    
    default:
        break;
    }
}

template<class Type, class Index>
mesh_2d<Type, Index>::mesh_2d(const std::string& path)
{
    read_su2(path);
}

template<class Type, class Index>
void mesh_2d<Type, Index>::read_su2(const std::string& path)
{
    std::ifstream mesh_file(path);

    if(mesh_file.is_open())
    {
        std::string pass;
        size_t count = 0;

        mesh_file >> pass >> pass >> pass >> count;

        elements_2d_type.resize(count);
        for(size_t i = 0; i < count; ++i)
        {
            int temp = 0;
            mesh_file >> temp;
            if(temp == 9)
            {
                elements_2d_type[i] = 3;
                elements.push_back(std::array<Index, 4>({}));
                mesh_file >> elements(i, 0)
                          >> elements(i, 3)
                          >> elements(i, 2)
                          >> elements(i, 1);
            }
            else if(temp == 23)
            {
                elements_2d_type[i] = 4;
                elements.push_back(std::array<Index, 8>({}));
                mesh_file >> elements(i, 0)
                          >> elements(i, 6)
                          >> elements(i, 4)
                          >> elements(i, 2)
                          >> elements(i, 7)
                          >> elements(i, 5)
                          >> elements(i, 3)
                          >> elements(i, 1);
            }
            mesh_file >> pass;
        }

        mesh_file >> pass >> count;
        nodes.resize(count, 2);
        for(size_t i = 0; i < count; ++i)
            mesh_file >> nodes(i, 0) >> nodes(i, 1) >> pass;

        mesh_file >> pass >> count;
        bounds.resize(count);
        elements_1d_type.resize(count);
        for(size_t b = 0; b < bounds.size(); ++b)
        {
            mesh_file >> pass >> pass;
            if(pass != std::string("Group_Of_All_Edges"))
            {
                mesh_file >> pass >> count;
                elements_1d_type[b].resize(count);
                for(size_t i = 0; i < count; ++i)
                {
                    int temp = 0;
                    mesh_file >> temp;
                    elements_1d_type[b][i] = temp - '0';
                    if(temp == 3)
                    {
                        elements_1d_type[b][i] = 0;
                        bounds[b].push_back(std::array<Index, 2>({}));
                        mesh_file >> bounds[b](i, 0) >> bounds[b](i, 1);
                    }

                    if(temp == 21)
                    {
                        elements_1d_type[b][i] = 1;
                        bounds[b].push_back(std::array<Index, 3>({}));
                        mesh_file >> bounds[b](i, 0) >> bounds[b](i, 1) >> bounds[b](i, 2);
                    }
                }
            }
            else
            {
                bounds.resize(bounds.size() - 1);
                elements_1d_type.resize(elements_1d_type.size() - 1);
            }
            
        }
        mesh_file.close();
    }
}

template<class Type, class Index>
void mesh_2d<Type, Index>::clear()
{
    elements.clear();
    nodes.clear();
    bounds.clear();
    neighbors.clear();
    elements_1d_type.clear();
    elements_2d_type.clear();
}

// Вычисление ориентированной площади треугольника на точках a, b, c.
template<class T>
T oriented_triangle_area(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c) noexcept {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Проверка отрезков AB и CD на предмет пересечения.
template<class T>
bool intersect(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c, const std::array<T, 2>& d) noexcept {
    static constexpr auto check_segment = [](T a, T b, T c, T d) {
        if (a > b) std::swap(a, b);
        if (c > d) std::swap(c, d);
        return std::max(a, c) <= std::min(b, d);
    };
    
    return check_segment(a[0], b[0], c[0], d[0]) && 
           check_segment(a[1], b[1], c[1], d[1]) &&
           oriented_triangle_area(a, b, c) * oriented_triangle_area(a, b, d) <= 0 &&
           oriented_triangle_area(c, d, a) * oriented_triangle_area(c, d, b) <= 0;
}

template<class Type, class Index>
bool mesh_2d<Type, Index>::check_intersect_with_bound(const std::array<Type, 2>& A, const std::array<Type, 2>& B)
{
    for(size_t b = 0; b < bounds.size(); ++b)
        for(size_t i = 0; i < bounds[b].rows(); ++i)
            for(size_t j = 0; j < bounds[b].cols(i)-1; ++j)
            {
                std::array<Type, 2> C = { coord(bounds[b](i, j  ), 0), coord(bounds[b](i, j  ), 1) },
                                    D = { coord(bounds[b](i, j+1), 0), coord(bounds[b](i, j+1), 1) };
                if(intersect(A, B, C, D))
                    return true;
            }
    return false;
}

template<class Type, class Index>
void mesh_2d<Type, Index>::find_neighbors_for_elements(const Type r)
{
    matrix<Type> centers(elements_count(), 2, 0.0); // Координаты центров элементов
    for(size_t el = 0; el < elements_count(); ++el)
    {
        const metamath::finite_element::element_2d_integrate_base<Type> *const e = element_2d(element_type(el));
        const Type x0 = el < 3 ? 1.0/3.0 : 0.0;
        for(size_t i = 0; i < 2; ++i)
            for(size_t j = 0; j < e->nodes_count(); ++j)
                centers(el, i) += coord(elements(el, j), i) * e->N(j, {x0, x0});
    }

    neighbors.resize(elements_count());
    const auto distance = [](const Type x1, const Type x2, const Type y1, const Type y2)
                          { return sqrt(metamath::power<2>(x1-x2) + metamath::power<2>(y1-y2)); };
    for(size_t elL = 0; elL < elements_count(); ++elL)
    {
        neighbors[elL].resize(0);
        neighbors[elL].reserve(elements_count());
        for(size_t elN = 0; elN < elements_count(); ++elN)
            if(distance(centers(elL, 0), centers(elN, 0), centers(elL, 1), centers(elN, 1)) < r &&
               !check_intersect_with_bound({centers(elL, 0), centers(elL, 1)}, {centers(elN, 0), centers(elN, 1)}))
                neighbors[elL].push_back(elN);
        neighbors[elL].shrink_to_fit();
    }
}

template<class Type, class Index>
void mesh_2d<Type, Index>::find_neighbors_for_nodes(const Type r)
{
    matrix<Type> centers(elements_count(), 2, 0.0); // Координаты центров элементов
    for(size_t el = 0; el < elements_count(); ++el)
    {
        const metamath::finite_element::element_2d_integrate_base<Type> *const e = element_2d(element_type(el));
        const Type x0 = el < 3 ? 1.0/3.0 : 0.0;
        for(size_t i = 0; i < 2; ++i)
            for(size_t j = 0; j < e->nodes_count(); ++j)
                centers(el, i) += coord(elements(el, j), i) * e->N(j, {x0, x0});
    }

    neighbors.resize(nodes_count());
    const auto distance = [](const Type x1, const Type x2, const Type y1, const Type y2)
                          { return sqrt(metamath::power<2>(x1-x2) + metamath::power<2>(y1-y2)); };
    for(size_t node = 0; node < nodes_count(); ++node)
    {
        neighbors[node].resize(0);
        neighbors[node].reserve(elements_count());
        for(size_t elN = 0; elN < elements_count(); ++elN)
            if(distance(coord(node, 0), centers(elN, 0), coord(node, 1), centers(elN, 1)) < r &&
               !check_intersect_with_bound({centers(node, 0), centers(node, 1)}, {centers(elN, 0), centers(elN, 1)}))
                neighbors[node].push_back(elN);
        neighbors[node].shrink_to_fit();
    }
}

template<class Type, class Index>
size_t mesh_2d<Type, Index>::nodes_count() const noexcept { return nodes.rows(); }
template<class Type, class Index>
size_t mesh_2d<Type, Index>::elements_count() const noexcept { return elements.rows(); }
template<class Type, class Index>
size_t mesh_2d<Type, Index>::boundary_groups_count() const noexcept { return bounds.size(); }

template<class Type, class Index>
Index mesh_2d<Type, Index>::node_number(const size_t el, const size_t i) const noexcept { return elements(el, i); }
template<class Type, class Index>
Type mesh_2d<Type, Index>::coord(const size_t node, const size_t component) const noexcept { return nodes(node, component); }

template<class Type, class Index>
const rows_different_sizes<Index>& mesh_2d<Type, Index>::boundary(const size_t b) const noexcept { return bounds[b]; }

template<class Type, class Index>
const std::vector<Index>& mesh_2d<Type, Index>::neighbor(const size_t el) const noexcept { return neighbors[el]; }

template<class Type, class Index>
uint8_t mesh_2d<Type, Index>::element_type(const size_t el) const noexcept { return elements_2d_type[el]; }
template<class Type, class Index>
const std::vector<uint8_t>& mesh_2d<Type, Index>::elements_on_bound_types(const size_t b) const noexcept { return elements_1d_type[b]; }

template<class Type, class Index>
metamath::finite_element::element_2d_integrate_base<Type>* mesh_2d<Type, Index>::element_2d(const size_t el) noexcept { return finite_elements_2d[el]; }
template<class Type, class Index>
metamath::finite_element::element_1d_integrate_base<Type>* mesh_2d<Type, Index>::element_1d(const size_t el) noexcept { return finite_elements_1d[el]; }

template<class Type, class Index>
const metamath::finite_element::element_2d_integrate_base<Type>* mesh_2d<Type, Index>::element_2d(const size_t el) const noexcept { return finite_elements_2d[el]; }
template<class Type, class Index>
const metamath::finite_element::element_1d_integrate_base<Type>* mesh_2d<Type, Index>::element_1d(const size_t el) const noexcept { return finite_elements_1d[el]; }

template<class Type, class Index>
void mesh_2d<Type, Index>::print_to_file(const std::string& path, const Type *x) const
{
    std::ofstream fout(path);
    fout.precision(20);
    for(size_t i = 0; i < nodes_count(); ++i)
        fout << coord(i, 0) << "," << coord(i, 1) << "," << x[i] << std::endl;
}

template<class Type, class Index>
mesh_2d<Type, Index>::~mesh_2d()
{
    for(size_t i = 0; i < finite_elements_1d.size(); ++i)
        delete finite_elements_1d[i];

    for(size_t i = 0; i < finite_elements_2d.size(); ++i)
        delete finite_elements_2d[i];
}

#endif