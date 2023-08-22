#include "determine_problem.hpp"

int main(const int argc, const char *const *const argv) {
    if (argc != 2) {
        std::cerr << "ERROR: Input format: [program name] path/to/config.json" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        using T = double;
        using I = int64_t;
        nonlocal::determine_problem<T, I>(nonlocal::config::parse_json(std::filesystem::path{argv[1]}));
    } catch (const std::exception& e) {
        using namespace std::literals;
        std::cerr << "ERROR: "s + e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "ERROR: unknown" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}