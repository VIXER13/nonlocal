#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <ostream>

namespace logger {

enum class log_level : uint8_t {
    ERROR,
    WARNING,
    INFO,
    DEBUG,
    TRACE,
    COUNT
};

struct stream_base {
    std::ostream& out;
    explicit stream_base(std::ostream& out) noexcept;
    virtual ~stream_base() noexcept = default;
};

class logger {
    const std::chrono::time_point<std::chrono::system_clock> _initial_time = std::chrono::system_clock::now();
    const log_level _level;
    std::unique_ptr<stream_base> _out;

    explicit logger(const log_level level, std::unique_ptr<stream_base>&& out);

public:
    log_level level() const noexcept;
    std::ostream& log(const log_level level = log_level::INFO);

    friend logger& get(const log_level level, std::unique_ptr<stream_base>&& init);
};

logger& get(const log_level level = log_level::INFO, std::unique_ptr<stream_base>&& init = nullptr);

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

#endif