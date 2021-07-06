#include <iostream>
#include <ostream>
#include "heat_equation_solver_1d.hpp"

namespace {

enum class element_type : uint8_t {
    LINEAR = 1,
    QUADRATIC = 2,
    QUBIC = 3
};

template<class T>
using finite_element_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

template<class T>
static finite_element_1d_ptr<T> make_element(const element_type type) {
    switch(type) {
        case element_type::LINEAR: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss1>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::linear>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUADRATIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss2>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::quadratic>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUBIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss3>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::qubic>;
            return std::make_unique<element_1d>(quadrature{});
        }

        default:
            throw std::logic_error{"Unknown element type " + std::to_string(int(type))};
    }
}

template<class Vector>
void save_as_csv(const std::string& path, const Vector& x) {
    std::ofstream csv{path};
    //csv.precision(std::numeric_limits<T>::max_digits10);
    const double h = 1. / (x.size() - 1);
    for(size_t i = 0; i < x.size(); ++i)
        csv << i * h << ',' << x[i] << '\n';
}

template<class T, uintmax_t p, uintmax_t q>
class polynomial final {
    T _r, _norm;

public:
    explicit polynomial(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = T{p} / (2 * _r * std::beta(T{1} / T{p}, T{q + 1}));
    }

    T radius() const noexcept { return _r; }

    T operator()(const T x, const T y) const noexcept {
        const T h = std::abs(x - y);
        using metamath::function::power;
        return h < _r ? _norm * power<q>(T{1} - power<p>(h / _r)) : T{0};
    }
};

template<class Func>
void test(const Func& func = [](){}) {
    func();
}

}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "run format: program_name <element_type> <elements_count> <section>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<mesh::mesh_1d<double, int>>(
            make_element<double>(element_type(std::stoi(argv[1]))),
            std::stoull(argv[2]), std::array{std::stod(argv[3]), std::stod(argv[4])}
        );

        nonlocal::solver_parameters<double> sol_parameters;
        sol_parameters.save_path = "./impulseLoc/";
        sol_parameters.time_interval[0] = 0;
        sol_parameters.time_interval[1] = 10;
        sol_parameters.steps = 10000;
        sol_parameters.save_freq = 1;

        nonlocal::heat::equation_parameters<double> parameters;
        parameters.p1 = 0.3;
        parameters.r = 0.1;
        mesh->calc_neighbours_count(parameters.r);
        nonlocal::heat::heat_equation_solver_1d<double, int> solver{mesh};

        //nonlocal::finite_element_solver_base_1d<double, int> solver{mesh};

//        solver.nonstationary(sol_parameters, parameters,
//            {
//                std::pair{
//                    nonlocal::boundary_condition_t::SECOND_KIND,
//                    [](const double t) { return 4 * t * t * std::exp(-2 * t); }
//                },
//                std::pair{
//                    nonlocal::boundary_condition_t::SECOND_KIND,
//                    [](const double t) { return 0; }
//                }
//            },
//            [](const double x) { return 0; },
//            [](const double x) { return 0; },
//            polynomial<double, 2, 1>{parameters.r}
//        );

        auto solution = solver.stationary(
            parameters,
            {
                std::pair{nonlocal::boundary_condition_t::FIRST_KIND, 1},
                std::pair{nonlocal::boundary_condition_t::FIRST_KIND, 0.5},
            },
            [](const double x) { return 0; },
            polynomial<double, 2, 1>{parameters.r}
        );
        save_as_csv("test.csv", solution);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}