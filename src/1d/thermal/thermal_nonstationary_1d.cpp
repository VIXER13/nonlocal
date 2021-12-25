#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver.hpp"
#include "influence_functions_1d.hpp"
#include <iostream>

int main(const int argc, const char *const *const argv) {
    if (argc < 6) {
        std::cerr << "run format: program_name <element_type> <elements_count> <p1> <r> <save_name>" << std::endl;
        return EXIT_FAILURE;
    }

    try {

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}