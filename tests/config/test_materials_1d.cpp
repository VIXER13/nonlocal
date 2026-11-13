#include <config/read_mechanical_parameters.hpp>
#include <config/read_mesh.hpp>
#include <config/read_thermal_parameters.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal::config;

const suite<"config_materials_1d"> _ = [] {
    "mesh_and_solver_use_the_same_model"_test = [] {
        const auto config = nlohmann::json::parse(R"({
            "materials": [
                {
                    "length": 1, "elements_count": 4,
                    "physical": { "conductivity": 1, "youngs_modulus": 1 },
                    "model": { "local_weight": 0.75, "nonlocal_radius": 0.5 },
                    "thermal_model": { "local_weight": 0.5, "nonlocal_radius": [0.75] },
                    "mechanical_model": { "local_weight": 0.25, "nonlocal_radius": 0.5, "search_radius": 0.25 }
                },
                {
                    "length": 1, "elements_count": 4,
                    "physical": { "conductivity": 2, "youngs_modulus": 2 },
                    "model": { "local_weight": 0.75, "nonlocal_radius": 0.5 }
                }
            ]
        })");
        const auto& materials = config["materials"];
        const auto thermal = read_thermal_parameters_1d<T>(materials, "materials");
        const auto mechanical = read_mechanical_parameters_1d<T>(materials, "materials");
        const auto thermal_mesh = read_mesh_1d<T>(config, {}, "thermal");
        const auto mechanical_mesh = read_mesh_1d<T>(config, {}, "mechanical");

        // The specific model wins. The search radius overrides the influence radius when given.
        expect(eq(thermal[0].model.local_weight, T{ 0.5 }));
        expect(eq(thermal_mesh->search_radius(0), T{ 0.75 }));
        expect(eq(mechanical[0].model.local_weight, T{ 0.25 }));
        expect(eq(mechanical_mesh->search_radius(0), T{ 0.25 }));

        // Both readers use the common model on the second segment.
        expect(eq(thermal[1].model.local_weight, T{ 0.75 }));
        expect(eq(thermal_mesh->search_radius(1), T{ 0.5 }));
        expect(eq(mechanical[1].model.local_weight, T{ 0.75 }));
        expect(eq(mechanical_mesh->search_radius(1), T{ 0.5 }));
    };
};

} // namespace
