#include "json_parser.hpp"

element_type get_element_type(const std::string& str) {
	if (str == "linear")
        return element_type::LINEAR;
    if (str == "quadratic")
		return element_type::QUADRATIC;
	if (str == "qubic")
		return element_type::QUBIC;
	if (str == "quartic")
		return element_type::QUARTIC;
	if (str == "quintic")
		return element_type::QUINTIC;
    throw std::domain_error("Unknown element type: " + str);    
}

boundary_kind_type get_boundary_type(const std::string& str) {
	if (str == "temperature")
		return boundary_kind_type::TEMPERATURE;
    if (str == "flux")
        return boundary_kind_type::FLUX;
	throw std::domain_error("Unknown boundary condition type: " + str); 
}

void read_json_not_template(not_template_parameters& par, const Json::Value& Jval) {
    par.path_to_save_temperature = Jval["save"]["folder"].asString() + "/" + Jval["save"]["temperature"].asString();
	par.path_to_save_flux = Jval["save"]["folder"].asString() + "/" + Jval["save"]["flux"].asString();
	par.path_to_save_info = Jval["save"]["folder"].asString() + "/" + Jval["save"]["info"].asString();
	par.steps_count = Jval["steps_count"].asInt();
	par.save_frequent = Jval["save_frequent"].asInt();

	par.element_order = Jval["element_order"].isInt() ? element_type(Jval["element_order"].asInt()) : 
	                    get_element_type(Jval["element_order"].asString());

	par.left_boundary_kind = Jval["boundaries"]["left"]["kind"].isInt() ? boundary_kind_type(Jval["boundaries"]["left"]["kind"].asInt()) : 
	                         get_boundary_type(Jval["boundaries"]["left"]["kind"].asString());

    par.right_boundary_kind = Jval["boundaries"]["right"]["kind"].isInt() ? boundary_kind_type(Jval["boundaries"]["right"]["kind"].asInt()) : 
	                          get_boundary_type(Jval["boundaries"]["right"]["kind"].asString());
}

void read_json_template(json_data& res, const Json::Value& Jval) {
    for (const auto& current_material : Jval["materials"]) {
        res.time_step = Jval["time_step"].asDouble();
        res.left_boundary_value = Jval["boundaries"]["left"]["value"].asDouble();
	    res.right_boundary_value = Jval["boundaries"]["right"]["value"].asDouble();
	    res.right_part = Jval["right_part"].asDouble();
        
		material buf;
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

Json::Value read_json_file(const std::string& json_file_name) {
	Json::Reader json_read;
	Json::Value json_value;
	std::ifstream json_file(json_file_name);
	json_read.parse(json_file, json_value);
	json_file.close();

	return json_value;
}

void get_data_from_json(json_data& res, const std::string& json_file_name) {
	Json::Value json_value = read_json_file(json_file_name);
    read_json_not_template(res.not_template, json_value);
    read_json_template(res, json_value);
}