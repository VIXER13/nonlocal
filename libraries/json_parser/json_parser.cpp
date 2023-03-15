#include "json_parser.hpp"

element_type Element_Type_String_To_Enum(const std::string& str) {
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
    throw std::domain_error("Unknown element type!: " + str);    
}

boundary_kind_type Boundary_Type_String_To_Enum(const std::string& str) {
	if (str == "temperature")
		return boundary_kind_type::TEMPERATURE;
    if (str == "flux")
        return boundary_kind_type::FLUX;
	throw std::domain_error("Unknown boundary condition type!: " + str); 
}

void read_json_not_template(not_template_parameters& par, const Json::Value& Jval) {
    par.path_to_save_temperature = Jval["save"]["folder"].asString() + "/" + Jval["save"]["temperature"].asString();
	par.path_to_save_flux = Jval["save"]["folder"].asString() + "/" + Jval["save"]["flux"].asString();
	par.path_to_save_info = Jval["save"]["folder"].asString() + "/" + Jval["save"]["info"].asString();
	par.steps_count = Jval["steps_count"].asInt();
	par.save_frequent = Jval["save_frequent"].asInt();

	if (Jval["element_order"].isInt())
		par.element_order = (element_type)Jval["element_order"].asInt();
	else
		par.element_order = Element_Type_String_To_Enum(Jval["element_order"].asString());

	if (Jval["boundaries"]["left"]["kind"].isInt())
		par.left_boundary_kind = (boundary_kind_type)Jval["boundaries"]["left"]["kind"].asInt();
	else
		par.left_boundary_kind = Boundary_Type_String_To_Enum(Jval["boundaries"]["left"]["kind"].asString());

	if (Jval["boundaries"]["right"]["kind"].isInt())
		par.right_boundary_kind = (boundary_kind_type)Jval["boundaries"]["right"]["kind"].asInt();
	else
		par.right_boundary_kind = Boundary_Type_String_To_Enum(Jval["boundaries"]["right"]["kind"].asString());
}
