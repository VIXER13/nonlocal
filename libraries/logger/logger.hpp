#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <ostream>

namespace logger {

enum class level : uint8_t {
    Off,
    Error,
    Warning,
    Info,
    Debug,
    Trace,
    Count
};

struct stream_base {
    std::ostream& out;
    explicit stream_base(std::ostream& out) noexcept;
    virtual ~stream_base() noexcept = default;
};

class logger {
    const std::chrono::time_point<std::chrono::system_clock> _initial_time = std::chrono::system_clock::now();
    const level _level;
    std::unique_ptr<stream_base> _out;

    explicit logger(const level level, std::unique_ptr<stream_base>&& out);

public:
    level log_level() const noexcept;
    std::chrono::time_point<std::chrono::system_clock> initial_time() const;

    friend logger& get(const level level, std::unique_ptr<stream_base>&& init);
    friend logger& error(std::unique_ptr<stream_base>&& init);
    friend logger& warning(std::unique_ptr<stream_base>&& init);
    friend logger& info(std::unique_ptr<stream_base>&& init);
    friend logger& debug(std::unique_ptr<stream_base>&& init);
    friend logger& trace(std::unique_ptr<stream_base>&& init);

    template<class T>
    friend logger& operator<<(logger& log, const T& value);
    friend logger& operator<<(logger& log, std::ostream&(*f)(std::ostream&));
};

logger& get(const level level = level::Info, std::unique_ptr<stream_base>&& init = nullptr);
logger& error(std::unique_ptr<stream_base>&& init = nullptr);
logger& warning(std::unique_ptr<stream_base>&& init = nullptr);
logger& info(std::unique_ptr<stream_base>&& init = nullptr);
logger& debug(std::unique_ptr<stream_base>&& init = nullptr);
logger& trace(std::unique_ptr<stream_base>&& init = nullptr);

template<class T>
logger& operator<<(logger& log, const T& value) {
    log._out->out << value;
    return log;
}

struct cout_stream : public stream_base {
    explicit cout_stream() noexcept;
};

struct cerr_stream : public stream_base {
    explicit cerr_stream() noexcept;
};

struct file_stream : public stream_base {
    std::fstream file;
    explicit file_stream(const std::filesystem::path& path);
};

}