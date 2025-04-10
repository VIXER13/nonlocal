#include "read_mechanical_parameters.hpp"

namespace nonlocal::config {

std::optional<std::bitset<2>> _mechanical_parameters_2d::read_null(const nlohmann::json& config, const std::string& path) {
    if (config.is_number())
        return std::nullopt;
    if (config.is_array() && config.size() == 2) {
        std::bitset<2> null;
        null[0] = config[0].is_null();
        null[1] = config[1].is_null();
        return null;
    }
    throw std::domain_error{"The mechanical parameter \"" + path + "\" "
                            "must be either a number in isotropic case, "
                            "or an array of size 2 in the orthotropic case."};
}

void _mechanical_parameters_2d::check_null_combinations(const std::optional<std::bitset<2>> is_null_youngs_modulus,
                                                        const std::optional<std::bitset<2>> is_null_poissons_ratio,
                                                        const std::string& path) {
    if (is_null_youngs_modulus && !is_null_poissons_ratio || !is_null_youngs_modulus && is_null_poissons_ratio)
        throw std::domain_error{"youngs_modulus or poissons_ratio is incorrectly specified for material \"" + path + "\". "
                                "For the isotropic case, both materials must be specified as a number, for the orthotropic case, as arrays of size 2."};
    if (is_null_youngs_modulus && is_null_poissons_ratio && is_null_youngs_modulus->count() + is_null_poissons_ratio->count() != 1)
        throw std::domain_error{"Error in the combination of parameters of youngs_modulus and poissons_ratio for material \"" + path + "\"."
                                "Three parameters must be defined, the fourth must be specified as null, since its calculation will be automatic."};
}

}