#include "make_element.hpp"
#include "mesh_1d.hpp"

#include <iostream>

namespace {
using T = double;
}

int main(const int argc, const char *const *const argv) {
    try {
        const nonlocal::mesh::mesh_1d<T> mesh{
            nonlocal::make_element<T>(nonlocal::element_type::QUINTIC),
            {
                {.length = 5.,  .elements = 100},
                {.length = 10., .elements = 50},
                {.length = 5.,  .elements = 200},
            }
        };

        std::cout << "sections count = " << mesh.sections_count() << std::endl;
        std::cout << "elements count = " << mesh.elements_count() << std::endl;
        std::cout << "nodes count = " << mesh.nodes_count() << std::endl;
        std::cout << "length = " << mesh.length() << std::endl;
        std::cout << std::endl;

        for(const size_t section : std::ranges::iota_view{size_t{0}, mesh.sections_count()}) {
            std::cout << "section " << section << std::endl;
            std::cout << "elements count = " << mesh.elements_count(section) << std::endl;
            std::cout << "nodes count = " << mesh.nodes_count(section) << std::endl;
            std::cout << "length = " << mesh.length(section) << std::endl;
            std::cout << "step = " << mesh.step(section) << std::endl;
            std::cout << "jacobian = " << mesh.jacobian(section) << std::endl;
            const std::array<T, 2> bounds = mesh.bounds(section);
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