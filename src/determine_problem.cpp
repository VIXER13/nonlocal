#include "determine_problem.hpp"

namespace nonlocal {

std::vector<std::string> _determine_problem::get_required_fields(const bool is_time_dependent) {
    if (is_time_dependent)
        return {"boundaries", "materials", "time"};
    else
        return {"boundaries", "materials"};
}
    
}