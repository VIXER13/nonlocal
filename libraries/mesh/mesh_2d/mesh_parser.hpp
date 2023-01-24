#ifndef NONLOCAL_MESH_PARSER_HPP
#define NONLOCAL_MESH_PARSER_HPP

#include <cstdint>

namespace nonlocal::mesh {

enum class mesh_format : uint8_t {
    SU2
};

template<class T, class I, mesh_format Format>
class mesh_parser;

}

#endif