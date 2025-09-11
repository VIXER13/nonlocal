#include "read_model.hpp"

#include <logger/logger.hpp>

namespace nonlocal::config {

std::string get_model_field(const nlohmann::json& config, const std::string& path_with_access, const std::string& prefix) {
    if (!prefix.empty())
        if (const std::string field = prefix + "_model"; config.contains(field))
            return field;
    if (config.contains("model"))
        return "model";
    logger::debug() << "Optional field \"" + path_with_access + "model\" is not specified." << std::endl;
    return "";
}

}