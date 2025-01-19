#include "determine_problem.hpp"

int main(int argc, char** argv) {
    if (argc != 2) {
        logger::error() << "Input format: [program name] path/to/config.json" << std::endl;
        return EXIT_FAILURE;
    }

    int result = EXIT_SUCCESS;
    try {
#ifdef MPI_BUILD
        MPI_Init(&argc, &argv);
#endif
        using T = double;
        using I = int64_t;
        logger::info() << "NonLocFEM started." << std::endl;
        nonlocal::determine_problem<T, I>(nonlocal::config::parse_json(std::filesystem::path{argv[1]}));
        logger::info() << "NonLocFEM finished." << std::endl;
    } catch (const std::exception& e) {
        logger::error() << e.what() << std::endl;
        result = EXIT_FAILURE;
    } catch (...) {
        logger::error() << "Unknown error." << std::endl;
        result = EXIT_FAILURE;
    }

#ifdef MPI_BUILD
    MPI_Finalize();
#endif
    return result;
}