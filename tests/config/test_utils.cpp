#include "test_utils_json.h"

#include <config/config_utils.hpp>
#include <config/read_model.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::config;

const suite<"config_utils"> _ = [] {
    using namespace std::literals;
    "append_access_sign"_test = []{
        expect(eq(append_access_sign(""), ""s));
        expect(eq(append_access_sign("", 0), "[0]"s));
        expect(eq(append_access_sign("", 10), "[10]"s));
        expect(eq(append_access_sign("test"), "test."s));
        expect(eq(append_access_sign("test", 0), "test[0]"s));
    };

    const nlohmann::json config = nlohmann::json::parse(test_utils_json_data);

    "check_required_fields"_test = [&config]{
        expect(throws<std::domain_error>([&config]{ check_required_fields(config["test_required_fields"], { "third" }); }));
        expect(nothrow([&config]{ check_required_fields(config["test_required_fields"], { "first", "second" }); }));
        expect(nothrow([&config]{ check_optional_fields(config["test_required_fields"], { "third" }); }));
    };

    "get_model_field"_test = [&config]{
        expect(eq(get_model_field(config["model_with_prefix"], "", "some_prefix"), "some_prefix_model"s));
        expect(eq(get_model_field(config["model_without_prefix"], "", "some_prefix"), "model"s));
        expect(eq(get_model_field(config["model_with_and_without_prefix"], "", ""), "model"s));
        expect(eq(get_model_field(config["no_model"], "", "some_prefix"), ""s));
        expect(eq(get_model_field(config["no_model"], "", ""), ""s));
    };
};

}