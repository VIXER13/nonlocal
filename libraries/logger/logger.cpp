#include "logger.hpp"

#include <parallel/MPI_utils.hpp>

#include <iostream>
#include <array>

namespace {

struct dummy_stream : public logger::stream_base {
    std::fstream empty;
    explicit dummy_stream() :
        logger::stream_base{empty} {}
};

}

namespace logger {

stream_base::stream_base(std::ostream& out) noexcept
    : out{out} {}

logger::logger(const level level, std::unique_ptr<stream_base>&& out)
    : _level{level}
    , _out{std::move(out)} {}

level logger::log_level() const noexcept {
    return _level;
}

std::chrono::time_point<std::chrono::system_clock> logger::initial_time() const {
    return _initial_time;
}

logger& get(const level level, std::unique_ptr<stream_base>&& init) {
    static constexpr auto levels = std::array{
        "OFF: ",
        "ERROR: ",
        "WARNING: ",
        "INFO: ",
        "DEBUG: ",
        "TRACE: "
    };
    static_assert(levels.size() == size_t(level::Count));

    static logger log{level, init ? std::move(init) : std::make_unique<cout_stream>()};
    if (uint8_t(level) <= uint8_t(log.log_level()) && level != level::Off) {
        const auto current_time = std::chrono::system_clock::now();
        const std::chrono::duration<double> duration = current_time - log._initial_time;
        log._out->out << '[' << duration.count() << "s] ";
        if (parallel::MPI_size() > 1)
            log._out->out << "PROCESS " << parallel::MPI_rank() << ' ';
        log._out->out << levels[size_t(level)];
        return log;
    }

    static logger dummy_log{level::Off, std::make_unique<dummy_stream>()};
    return dummy_log;
}

logger& error(std::unique_ptr<stream_base>&& init) {
    return get(level::Error, std::move(init));
}

logger& warning(std::unique_ptr<stream_base>&& init) {
    return get(level::Warning, std::move(init));
}

logger& info(std::unique_ptr<stream_base>&& init) {
    return get(level::Info, std::move(init));
}

logger& debug(std::unique_ptr<stream_base>&& init) {
    return get(level::Debug, std::move(init));
}

logger& trace(std::unique_ptr<stream_base>&& init) {
    return get(level::Trace, std::move(init));
}

logger& operator<<(logger& log, std::ostream&(*f)(std::ostream&)) {
    f(log._out->out);
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