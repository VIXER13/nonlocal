#include "logger.hpp"

#include <parallel/MPI_utils.hpp>

#include <algorithm>
#include <array>
#include <iostream>

namespace {
struct dummy_stream : public logger::stream_base {
    std::ostream empty = std::ostream{nullptr};
    explicit dummy_stream() :
        logger::stream_base{empty} {}
};

logger::level g_log_level = logger::level::Info;
}

namespace logger {

stream_base::stream_base(std::ostream& out) noexcept
    : out{out} {}

logger::logger(std::unique_ptr<stream_base>&& out)
    : _out{std::move(out)} {}

std::chrono::time_point<std::chrono::system_clock> logger::initial_time() const {
    return _initial_time;
}

void set_log_level(level l) noexcept {
    g_log_level = l;
}

std::optional<level> parse_level(std::string_view name) {
    // Case insensitive comparison
    const auto compare = [](std::string_view a, std::string_view b) {
        constexpr auto cmp = [](char c1, char c2) { return std::tolower(c1) == std::tolower(c2); };
        if (a.size() != b.size()) return false;
        return std::equal(a.begin(), a.end(), b.begin(), cmp);
    };

    if (compare(name, "off")) return level::Off;
    if (compare(name, "error")) return level::Error;
    if (compare(name, "warning")) return level::Warning;
    if (compare(name, "info")) return level::Info;
    if (compare(name, "debug")) return level::Debug;
    if (compare(name, "trace")) return level::Trace;
    return std::nullopt;
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

    static logger log{init ? std::move(init) : std::make_unique<cout_stream>()};
    if (uint8_t(level) <= uint8_t(g_log_level) && level != level::Off) {
        const auto current_time = std::chrono::system_clock::now();
        const std::chrono::duration<double> duration = current_time - log._initial_time;
        log._out->out << '[' << duration.count() << "s] ";
        if (parallel::MPI_size() > 1)
            log._out->out << "PROCESS " << parallel::MPI_rank() << ' ';
        log._out->out << levels[size_t(level)];
        return log;
    }

    static logger dummy_log{std::make_unique<dummy_stream>()};
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