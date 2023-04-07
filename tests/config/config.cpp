#include "test_config_json.h"

#include "thermal_config_data.hpp"

#include <gtest/gtest.h>

namespace {

const Json::Value tests = nonlocal::config::read_json(test_config_json_data, std::next(test_config_json_data, test_config_json_size));

}

TEST(configs, test_embedded)
{
    EXPECT_EQ(tests["test_1"]["x"].asInt(), 5);
}