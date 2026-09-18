#include "determine_problem.hpp"

#include <string_view>

void help(std::string_view program_name) {
    std::cout << "Usage: " << program_name << " [--log-level=<level>] <config_file_path>" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --log-level=<level>   Set the log level (off, error, warning, info, debug, trace)" << std::endl;
    std::cout << "  --help, -h            Show this help message" << std::endl;
}

int main(int argc, char** argv) {
    if(argc == 1 or argc > 3) {
        help(argv[0]);
        return EXIT_FAILURE;
    }

    int path_pos = 0;
    for(int i = 1; i < argc; ++i) {
        std::string_view arg{argv[i]};
        if (arg == "--help" || arg == "-h") {
            help(argv[0]);
            return EXIT_SUCCESS;
        } else if (arg.starts_with("--log-level=")) {
            std::string_view level_str = arg.substr(12);
            if (auto level_result = logger::parse_level(level_str); level_result) {
                logger::set_log_level(*level_result);
            } else {
                logger::error() << "Invalid log level: " << level_str << std::endl;
                return EXIT_FAILURE;
            }
        } else  {
            path_pos = i;
        }
    }

    if(path_pos == 0) {
        logger::error() << "No config file path provided." << std::endl;
        return EXIT_FAILURE;
    }

    int result = EXIT_SUCCESS;
    try {
#ifdef MPI_BUILD
        MPI_Init(&argc, &argv);
#endif
        using T = double;
        logger::info() << "NonLocFEM started." << std::endl;
        nonlocal::determine_problem<T>(nonlocal::config::parse_json(std::filesystem::path{argv[path_pos]}));
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