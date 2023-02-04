#ifndef NONLOCAL_SOLUTION_2D_HPP
#define NONLOCAL_SOLUTION_2D_HPP

#include "mesh_2d.hpp"
#include "../solvers_constants.hpp"
#include "../equation_parameters.hpp"

namespace nonlocal {

template<class T, class I>
class solution_2d {
    const std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    const std::unordered_map<std::string, model_parameters<2, T>> _models;

protected:
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    explicit solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                         const std::unordered_map<std::string, model_parameters<2, T>>& models);

    // template<class Callback>
    // void calc_nonlocal(const Callback& callback) const;

    void save_scalars(std::ofstream& output, const std::vector<T>& x, const std::string_view name) const;
    void save_vectors(std::ofstream& output, const std::array<std::vector<T>, 2>& vec, const std::string_view name) const;

public:
    virtual ~solution_2d() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    const model_parameters<2, T>& model(const std::string& group) const;

    virtual void save_as_vtk(std::ofstream& output) const;
};

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
solution_2d<T, I>::solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                               const std::unordered_map<std::string, model_parameters<2, T>>& models)
    : _mesh{mesh}
    , _models{models} {}

template<class T, class I>
const mesh::mesh_2d<T, I>& solution_2d<T, I>::mesh() const {
    return *_mesh;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_2d<T, I>>& solution_2d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
const model_parameters<2, T>& solution_2d<T, I>::model(const std::string& group) const {
    return _models.at(group);
}

// template<class T, class I>
// template<class Callback>
// void solution_2d<T, I>::calc_nonlocal(const Callback& callback) const {
//     const auto process_nodes = mesh().process_nodes();
//     std::vector<bool> neighbors(mesh().container().elements_count(), true);
// #pragma omp parallel for default(none) shared(process_nodes) firstprivate(neighbors, callback)
//     for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
//         std::fill(neighbors.begin(), neighbors.end(), true);
//         for(const I eL : mesh().elements(node))
//             for(const I eNL : mesh().neighbours(eL))
//                 if (neighbors[eNL]) {
//                     callback(eNL, node);
//                     neighbors[eNL] = false;
//                 }
//     }
// }

template<class T, class I>
void solution_2d<T, I>::save_scalars(std::ofstream& output, const std::vector<T>& x, const std::string_view name) const {
    output << "SCALARS " << name << ' ' << mesh::vtk_data_type<T> << " 1\n"
           << "LOOKUP_TABLE default\n";
    for(const T val : x)
        output << val << '\n';
}

template<class T, class I>
void solution_2d<T, I>::save_vectors(std::ofstream& output, const std::array<std::vector<T>, 2>& vec, const std::string_view name) const {
    output << "VECTORS " << name << ' ' << mesh::vtk_data_type<T> << '\n';
    for(const size_t i : std::ranges::iota_view{0u, mesh().container().nodes_count()})
        output << vec[X][i] << ' ' << vec[Y][i] << " 0\n";
}

template<class T, class I>
void solution_2d<T, I>::save_as_vtk(std::ofstream& output) const {
    output.precision(std::numeric_limits<T>::max_digits10);
    mesh::utils::save_as_vtk(output, mesh().container());
    output << "POINT_DATA " << mesh().container().nodes_count() << '\n';
}

}

#endif