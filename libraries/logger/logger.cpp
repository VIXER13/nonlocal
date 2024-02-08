#include "logger.hpp"

#include "MPI_utils.hpp"

#include <iostream>
#include <array>

namespace logger {

stream_base::stream_base(std::ostream& out) noexcept
    : out{out} {}

logger::logger(const log_level level, std::unique_ptr<stream_base>&& out)
    : _level{level}
    , _out{std::move(out)} {}

log_level logger::level() const noexcept {
    return _level;
}

std::ostream& logger::log(const log_level level) {
    static constexpr auto levels = std::array{
        "ERROR: ",
        "WARNING: ",
        "INFO: ",
        "DEBUG: ",
        "TRACE: "
    };
    static_assert(levels.size() == size_t(log_level::COUNT));
    if (uint8_t(level) <= uint8_t(_level)) {
        const auto current_time = std::chrono::system_clock::now();
        const std::chrono::duration<double> duration = current_time - _initial_time;
        _out->out << '[' << duration.count() << "s] ";
        if (parallel::MPI_size() > 1)
            _out->out << "PROCESS " << parallel::MPI_rank() << ' ';
        _out->out << levels[size_t(level)];
        return _out->out;
    }
    static std::fstream empty;
    return empty;
}

logger& get(const log_level level, std::unique_ptr<stream_base>&& init) {
    static logger log{level, init ? std::move(init) : std::make_unique<cout_stream>()};
    return log;
}

cout_stream::cout_stream() noexcept
    : stream_base{std::cout} {}

cerr_stream::cerr_stream() noexcept
    : stream_base{std::cerr} {}

file_stream::file_stream(const std::filesystem::path& path)
    : stream_base{file}
    , file{path, std::ofstream::out} {}

}