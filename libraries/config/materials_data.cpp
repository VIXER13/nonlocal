#include "materials_data.hpp"

namespace nonlocal::config {

[[noreturn]] std::string _materials_data::throw_error(const std::string& field, const std::string& type) {
    throw std::domain_error{
        field.empty() ? 
        "materials_data initialization requires the initializing config to be an nonempty " + type :
        "materials_data initialization requires field \"" + field + "\" to be an nonempty " + type
    };
}

}