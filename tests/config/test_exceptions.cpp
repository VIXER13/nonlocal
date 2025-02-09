#include "test_exceptions_json.h"

#include <config/read_mechanical_boundary_conditions.hpp>
#include <config/read_model.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::config;

template<class T>
void expect_throws(T&& reader) {
    expect(throws<std::domain_error>(reader));
}

template<class T>
void expect_nothrows(T&& reader) {
    expect(nothrow(reader));
}

const suite<"config_exceptions"> _ = [] {
    using T = double;
    const nlohmann::json config = nlohmann::json::parse(test_exceptions_json_data);

    expect_throws([&config]{read_model_1d<T>(config["model_1d_local_weight_is_missed"], "");});
    expect_throws([&config]{read_model_1d<T>(config["model_1d_nonlocal_radius_is_missed"], "");});
    expect_throws([&config]{read_model_1d<T>(config["model_1d_local_weight_is_negative"], "");});
    expect_throws([&config]{read_model_1d<T>(config["model_1d_local_weight_is_zero"], "");});
    expect_throws([&config]{read_model_1d<T>(config["model_1d_local_weight_greater_1"], "");});
    expect_throws([&config]{read_model_1d<T>(config["model_1d_nonlocal_radius_is_negative"], "");});
    expect_nothrows([&config]{read_model_1d<T>(config["model_1d_ok"], "");});

    expect_throws([&config]{read_model_2d<T>(config["model_2d_local_weight_is_missed"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_nonlocal_radius_is_missed"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_nonlocal_radius_wrong_dimension"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_local_weight_is_negative"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_local_weight_is_zero"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_local_weight_greater_1"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_nonlocal_radius_is_negative"], "");});
    expect_throws([&config]{read_model_2d<T>(config["model_2d_nonlocal_one_of_the_radius_is_negative"], "");});
    expect_nothrows([&config]{read_model_2d<T>(config["model_2d_ok"], "");});
    expect_nothrows([&config]{read_model_2d<T>(config["model_2d_two_radii_ok"], "");});

    expect_nothrows([&config]{read_mechanical_boundaries_conditions_2d<T>(config["mechanical_boundaries_conditions_2d_empty"], "");});
    expect_throws([&config]{read_mechanical_boundaries_conditions_2d<T>(config["mechanical_boundaries_conditions_2d_empty_condition"], "");});
    expect_throws([&config]{read_mechanical_boundaries_conditions_2d<T>(config["mechanical_boundaries_conditions_2d_wrong_dimension_1"], "");});
    expect_throws([&config]{read_mechanical_boundaries_conditions_2d<T>(config["mechanical_boundaries_conditions_2d_wrong_dimension_3"], "");});
    expect_nothrows([&config]{read_mechanical_boundaries_conditions_2d<T>(config["mechanical_boundaries_conditions_2d_ok"], "");});
};

}