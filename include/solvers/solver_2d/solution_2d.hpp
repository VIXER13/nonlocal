#ifndef NONLOCAL_SOLUTION_2D_HPP
#define NONLOCAL_SOLUTION_2D_HPP

#include "mesh_2d.hpp"
#include "../solvers_constants.hpp"

namespace nonlocal {

template<class T, class I>
class solution_2d {
protected:
    using Influence_Function = std::function<T(const std::array<T, 2>& x, const std::array<T, 2>& y)>;

private:
    static constexpr std::string_view vtk_data_type = std::is_same_v<T, float> ? "float" : "double";

    const std::shared_ptr<mesh::mesh_proxy<T, I>> _mesh_proxy;
    const T _local_weight = T{1};
    const Influence_Function _influence_function = [](const std::array<T, 2>&, const std::array<T, 2>&) constexpr noexcept { return T{}; };

protected:
    explicit solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    explicit solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                         const T local_weight, const Influence_Function& influence_function);

    template<class Callback>
    void calc_nonlocal(const Callback& callback) const;

    void save_scalars(std::ofstream& output, const std::vector<T>& x, const std::string_view name) const;
    void save_vectors(std::ofstream& output, const std::array<std::vector<T>, 2>& vec, const std::string_view name) const;

public:
    virtual ~solution_2d() noexcept = default;

    const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy() const noexcept;
    T local_weight() const noexcept;
    const Influence_Function& influence_function() const noexcept;

    virtual void save_as_vtk(std::ofstream& output) const;
};

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _mesh_proxy{mesh_proxy} {}

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                               const T local_weight, const Influence_Function& influence_function)
    : _mesh_proxy{mesh_proxy}
    , _local_weight{local_weight}
    , _influence_function{influence_function} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_proxy<T, I>>& solution_2d<T, I>::mesh_proxy() const noexcept {
    return _mesh_proxy;
}

template<class T, class I>
T solution_2d<T, I>::local_weight() const noexcept {
    return _local_weight;
}

template<class T, class I>
const typename solution_2d<T, I>::Influence_Function& solution_2d<T, I>::influence_function() const noexcept {
    return _influence_function;
}

template<class T, class I>
template<class Callback>
void solution_2d<T, I>::calc_nonlocal(const Callback& callback) const {
    std::vector<bool> neighbors(mesh_proxy()->mesh().elements_count(), true);
#pragma omp parallel for default(none) firstprivate(neighbors, callback)
    for(size_t node = mesh_proxy()->first_node(); node < mesh_proxy()->last_node(); ++node) {
        std::fill(neighbors.begin(), neighbors.end(), true);
        for(const I eL : mesh_proxy()->nodes_elements_map(node))
            for(const I eNL : mesh_proxy()->neighbors(eL))
                if (neighbors[eNL]) {
                    callback(eNL, node);
                    neighbors[eNL] = false;
                }
    }
}

template<class T, class I>
void solution_2d<T, I>::save_scalars(std::ofstream& output, const std::vector<T>& x, const std::string_view name) const {
    output << "SCALARS " << name << ' ' << vtk_data_type << " 1\n"
           << "LOOKUP_TABLE default\n";
    for(const T val : x)
        output << val << '\n';
}

template<class T, class I>
void solution_2d<T, I>::save_vectors(std::ofstream& output, const std::array<std::vector<T>, 2>& vec, const std::string_view name) const {
    output << "VECTORS " << name << ' ' << vtk_data_type << '\n';
    for(const size_t i : std::views::iota(size_t{0}, mesh_proxy()->mesh().nodes_count()))
        output << vec[X][i] << ' ' << vec[Y][i] << " 0\n";
}

template<class T, class I>
void solution_2d<T, I>::save_as_vtk(std::ofstream& output) const {
    output.precision(std::numeric_limits<T>::max_digits10);
    mesh_proxy()->mesh().save_as_vtk(output);
    output << "POINT_DATA " << mesh_proxy()->mesh().nodes_count() << '\n';
}

}

#endif