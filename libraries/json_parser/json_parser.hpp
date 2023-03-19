#include <fstream>
#include <string>
#include <vector>
#include <json\json.h>

using T = double;

enum class element_type : uint8_t {
	LINEAR = 1,
	QUADRATIC = 2,
	QUBIC = 3,
	QUARTIC = 4,
	QUINTIC = 5
};

enum class boundary_kind_type : uint8_t {
	TEMPERATURE = 1,
	FLUX = 2
};

struct material {
	T length = 1.0;
	size_t elements = 100;
	T local_weight = 1.0;
	T nonlocal_radius = 0.0;
	T search_radius = 0.0;
	T conductivity = 1.0;
	T capacity = 1.0;
	T density = 1.0;
};

struct not_template_parameters {
	std::string path_to_save_temperature;
	std::string path_to_save_flux;
    std::string path_to_save_info;
	uint64_t steps_count = 100;
	uint64_t save_frequent = 1;
	element_type element_order = element_type::LINEAR;
	boundary_kind_type left_boundary_kind = boundary_kind_type::FLUX;
	boundary_kind_type right_boundary_kind = boundary_kind_type::FLUX;
};

struct json_data {
	not_template_parameters not_template;
	T time_step = 0.01;
	T left_boundary_value = 0.0;
	T right_boundary_value = 0.0;
	T right_part = 0.0;
	std::vector<material> materials;
};

element_type get_element_type(const std::string& str);
boundary_kind_type get_boundary_type(const std::string& str);

void read_json_not_template(not_template_parameters& par, const Json::Value& Jval);
void read_json_template(json_data& res, const Json::Value& Jval);

Json::Value read_json_file(const std::string& json_file_name);

void get_data_from_json(json_data& res, const std::string& json_file_name);