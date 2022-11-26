#include "make_element.hpp"
#include "mesh_1d.hpp"
#include "thermal/thermal_conductivity_matrix_1d.hpp"
#include "thermal/heat_equation_parameters_1d.hpp"

#include <iostream>

namespace {
using T = double;
using I = int64_t;
}

int main(const int argc, const char *const *const argv) {
    try {
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
            nonlocal::make_element<T>(nonlocal::element_type::QUADRATIC),
            std::vector{
                nonlocal::mesh::segment_data{.length = 5., .elements = 5},
                nonlocal::mesh::segment_data{.length = 5., .elements = 5},
                nonlocal::mesh::segment_data{.length = 5., .elements = 5},
                //nonlocal::mesh::segment_data{.length = 10., .elements = 5},
                //nonlocal::mesh::segment_data{.length = 5.,  .elements = 20},
            }
        );

        std::cout << "segments count = " << mesh->segments_count() << std::endl;
        std::cout << "elements count = " << mesh->elements_count() << std::endl;
        std::cout << "nodes count = " << mesh->nodes_count() << std::endl;
        std::cout << "length = " << mesh->length() << std::endl;
        std::cout << std::endl;

        for(const size_t segment : std::ranges::iota_view{size_t{0}, mesh->segments_count()}) {
            std::cout << "segment " << segment << std::endl;
            std::cout << "elements count = " << mesh->elements_count(segment) << std::endl;
            std::cout << "nodes count = " << mesh->nodes_count(segment) << std::endl;
            std::cout << "length = " << mesh->length(segment) << std::endl;
            std::cout << "step = " << mesh->step(segment) << std::endl;
            std::cout << "jacobian = " << mesh->jacobian(segment) << std::endl;

            std::cout << "elements = ";
            for(const size_t e : mesh->segment_elements(segment))
                std::cout << e << ' ';
            std::cout << std::endl;

            std::cout << "nodes = ";
            for(const size_t node : mesh->segment_nodes(segment))
                std::cout << node << ' ';
            std::cout << std::endl;

            const std::array<T, 2> bounds = mesh->bounds(segment);
            std::cout << "bounds = [" << bounds.front() << ',' << bounds.back() << ']' << std::endl;
            std::cout << std::endl;
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}