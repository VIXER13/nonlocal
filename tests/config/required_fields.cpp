#include "tests_config_utils.hpp"
#include "required_fields_json.h"

#include "nonlocal_config.hpp"

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace unit_tests;
using namespace nonlocal::config;

template<class T>
auto expect_throw(const nlohmann::json& name) {
    return [&name] {
        expect(throws([&name] { T{name}; }));
    };
}

template<class T>
auto expect_no_throw(const nlohmann::json& name) {
    return [&name] {
        expect(nothrow([&name] { T{name}; }));
    };
}

const suite<"config_required_fields"> _ = [] {
    const nlohmann::json config = nlohmann::json::parse(required_fields_json_data);
};

}