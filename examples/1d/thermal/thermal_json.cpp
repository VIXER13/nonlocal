#include "json_parser.hpp"

int main(const int argc, const char *const *const argv)
{
    json_data<double> data;

    get_data_from_json(data, argv[1]);

    std::cout << data.element_order << std::endl;
}