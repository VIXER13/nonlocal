#include "determine_problem.hpp"

#include <argparse/cli.hpp>

#include <string_view>

int main(int argc, char** argv) {
    cli cmd;
    auto log_level = cmd.arg("--log-level", "Set the log level.").alias("-ll").cast(logger::parse_level).type<std::string>();
    auto cfg_path = cmd.pos("Path to the config file").type<std::string>();
    cmd.parse(argc, argv);

    if (!cfg_path.has_value()) {
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
        nonlocal::determine_problem<T>(nonlocal::config::parse_json(std::filesystem::path{ *cfg_path }));
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