#include <solvers/solver_1d/base/linear_solver_1d.hpp>
#include <tests/solvers_1d/utils/solution_checks.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using I = int64_t;
using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using sparse_matrix = Eigen::SparseMatrix<T, Eigen::RowMajor, I>;
using namespace nonlocal::solver_1d;
using namespace nonlocal::unit_tests;
using namespace boost::ut;

const suite<"linear_solver_1d"> _ = [] {
    "cholesky_fallback_and_singular_system"_test = [] {
        sparse_matrix matrix{ 2, 2 };
        matrix.insert(0, 1) = T{ 1 };
        vector right_part(2);
        right_part << T{ 3 }, T{ 2 };
        // The full symmetric matrix [0 1; 1 0] is indefinite: LLT fails, but LU gives (2, 3).
        static constexpr auto Expected = std::array{ T{ 2 }, T{ 3 } };
        check_values(solve_linear_system_1d(matrix, right_part, {}), Expected);
        matrix.setZero();
        expect(throws<std::domain_error>([&] { solve_linear_system_1d(matrix, right_part, {}); }));
    };
};

} // namespace
