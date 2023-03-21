#include <fstream>
#include <string>
#include <vector>
#include <json\json.h>

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
	double length = 1.0;
	size_t elements = 100;
	double local_weight = 1.0;
	double nonlocal_radius = 0.0;
	double search_radius = 0.0;
	double conductivity = 1.0;
	double capacity = 1.0;
	double density = 1.0;
};

struct json_data {
	std::string path_to_save_temperature;
	std::string path_to_save_flux;
    std::string path_to_save_info;
	uint64_t steps_count = 100;
	uint64_t save_frequent = 1;
	element_type element_order = element_type::LINEAR;
	boundary_kind_type left_boundary_kind = boundary_kind_type::FLUX;
	boundary_kind_type right_boundary_kind = boundary_kind_type::FLUX;
	double time_step = 0.01;
	double left_boundary_value = 0.0;
	double right_boundary_value = 0.0;
	double right_part = 0.0;
	std::vector<material> materials;
};

element_type get_element_type(const std::string& str);
boundary_kind_type get_boundary_type(const std::string& str);

void read_json_parameters(json_data& res, const Json::Value& Jval);
void read_json_materials(json_data& res, const Json::Value& Jval);

Json::Value read_json_file(const std::string& json_file_name);

void get_data_from_json(json_data& res, const std::string& json_file_name);