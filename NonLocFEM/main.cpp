#include "determine_problem.hpp"

int main(int argc, char** argv) {
    if (argc != 2) {
        logger::get().log(logger::log_level::ERROR) << "Input format: [program name] path/to/config.json" << std::endl;
        return EXIT_FAILURE;
    }

    int result = EXIT_SUCCESS;
    try {
#ifdef MPI_BUILD
        MPI_Init(&argc, &argv);
#endif
        using T = double;
        using I = int64_t;
        logger::get().log() << "NonLocFEM started." << std::endl;
        nonlocal::determine_problem<T, I>(nonlocal::config::parse_json(std::filesystem::path{argv[1]}));
        logger::get().log() << "NonLocFEM finished." << std::endl;
    } catch (const std::exception& e) {
        logger::get().log(logger::log_level::ERROR) << e.what() << std::endl;
        result = EXIT_FAILURE;
    } catch (...) {
        logger::get().log(logger::log_level::ERROR) << "Unknown error." << std::endl;
        result = EXIT_FAILURE;
    }

#ifdef MPI_BUILD
    MPI_Finalize();
#endif
    return result;
}