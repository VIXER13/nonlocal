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

template<typename T>
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
	boundary_kind_type left_boundary_kind = boundary_kind_type::TEMPERATURE;
	boundary_kind_type right_boundary_kind = boundary_kind_type::TEMPERATURE;
};

template<typename T>
struct json_data {
	not_template_parameters not_template;
	T time_step = 0.01;
	T left_boundary_value = 0.0;
	T right_boundary_value = 0.0;
	T right_part = 0.0;
	std::vector<material<T>> materials;
};

element_type Element_Type_String_To_Enum(const std::string& str);
boundary_kind_type Boundary_Type_String_To_Enum(const std::string& str);

void read_json_not_template(not_template_parameters& par, const Json::Value& Jval);

template<typename T>
void read_json_template(json_data<T>& res, const Json::Value& Jval) {
    for (const auto& current_material : Jval["materials"]) {
        res.time_step = Jval["time_step"].asDouble();
        res.left_boundary_value = Jval["boundaries"]["left"]["value"].asDouble();
	    res.right_boundary_value = Jval["boundaries"]["right"]["value"].asDouble();
	    res.right_part = Jval["right_part"].asDouble();
        
		material<T> buf;
		buf.length = current_material["length"].asDouble();
		buf.elements = current_material["elements"].asInt();
		buf.conductivity = current_material["parameters"]["physical"]["conductivity"].asDouble();
		buf.capacity = current_material["parameters"]["physical"]["capacity"].asDouble();
		buf.density = current_material["parameters"]["physical"]["density"].asDouble();

		if (current_material["parameters"].isMember("model")) {
			buf.local_weight = current_material["parameters"]["model"]["local_weight"].asDouble();
			buf.nonlocal_radius = current_material["parameters"]["model"]["nonlocal_radius"].asDouble();
			buf.search_radius = current_material["parameters"]["model"]["search_radius"].asDouble();
		}

		res.materials.push_back(buf);
	}
}

template<typename T>
void get_data_from_json(json_data<T>& res, const std::string& json_file_name) {
	Json::Reader json_read;
	Json::Value json_value;
	std::ifstream json_file(json_file_name);
	json_read.parse(json_file, json_value);
	json_file.close();

    read_json_not_template(res.not_template, json_value);
    read_json_template(res, json_value);
}