#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <ostream>
#include <fstream>
#include <filesystem>
#include <memory>

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
    const log_level _level;
    std::unique_ptr<stream_base> _out;

    explicit logger(const log_level level, std::unique_ptr<stream_base>&& out);

public:
    log_level level() const noexcept;
    std::ostream& log(const log_level level);

    friend logger& get(const log_level level, std::unique_ptr<stream_base>&& init);
};

logger& get(const log_level level = log_level::ERROR, std::unique_ptr<stream_base>&& init = nullptr);

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