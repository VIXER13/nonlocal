#include "create_matrix.hpp"

#include <solvers/slae/stable_biconjugate_gradient.hpp>

#include <solvers/slae/ilu0_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

void make_diagonal_dominant(sparse_matrix<T>& matrix) {
    for(size_t row = 0; row < matrix.rows(); ++row) {
        T sum = T{0};
        for(const size_t col : matrix.portrait.indices_range(row))
            if (col != row)
                sum += std::abs(matrix(row, col));
        if (T& val = matrix(row, row); val < sum)
            val = sum + static_cast<T>(std::rand()) / RAND_MAX;
    }
}

sparse_matrix<T> random_general_matrix(const size_t size) {
    sparse_matrix<T> matrix{size, size};
    for(size_t row = 0; row < size; ++row) {
        std::vector<bool> cols(size, false);
        for(size_t col = 0; col < size; ++col) {
            if (row == col || std::rand() % 2 == 0 && !cols[col]) {
                cols[col] = true;
                ++matrix.portrait.shifts[row + 1];
                matrix.portrait.indices.push_back(col);
                matrix.values.push_back(static_cast<T>(std::rand()) / RAND_MAX);
            }
        }
    }
    matrix.portrait.accumulate_shifts();
    make_diagonal_dominant(matrix);
    return matrix;
}

std::vector<T> random_vector(const size_t size) {
    std::vector<T> vec(size);
    for(auto& val : vec)
        val = static_cast<T>(std::rand()) / RAND_MAX;
    return vec;
}

void make_diagonal_dominant(sparse_matrix<square_matrix<T, 2>>& matrix) {
    for(size_t row = 0; row < matrix.rows(); ++row) {
        std::array<T, 2> sum{};
        for(const size_t col : matrix.portrait.indices_range(row))
            if (col != row) {
                sum[0] += std::abs(matrix(row, col)[0][0]) + std::abs(matrix(row, col)[0][1]);
                sum[1] += std::abs(matrix(row, col)[1][0]) + std::abs(matrix(row, col)[1][1]);
            }
        square_matrix<T, 2>& val = matrix(row, row);
        //if (val[0][0] < sum[0])
            val[0][0] = 1.001 * sum[0];
        //if (val[1][1] < sum[1])
            val[1][1] = 1.001 * sum[1];
    }
}

sparse_matrix<square_matrix<T, 2>> random_block_general_matrix(const size_t size) {
    sparse_matrix<square_matrix<T, 2>> matrix{size, size};
    constexpr auto rand = []() { return (static_cast<T>(std::rand()) - RAND_MAX/2); };
    for(size_t row = 0; row < size; ++row) {
        std::vector<bool> cols(size, false);
        for(size_t col = 0; col < size; ++col) {
            if (row == col || std::rand() % 2 == 0 && !cols[col]) {
                cols[col] = true;
                ++matrix.portrait.shifts[row + 1];
                matrix.portrait.indices.push_back(col);
                matrix.values.push_back(square_matrix<T, 2>{
                    rand(), rand(),
                    rand(), rand()
                });
            }
        }
    }
    matrix.portrait.accumulate_shifts();
    make_diagonal_dominant(matrix);
    return matrix;
}

std::vector<std::array<T, 2>> random_block_vector(const size_t size) {
    std::vector<std::array<T, 2>> vec(size);
    for(auto& val : vec) {
        val[0] = static_cast<T>(std::rand());
        val[1] = static_cast<T>(std::rand());
    }
    return vec;
}

suite<"stable_biconjugate_gradient"> _ = [] {
    "scalar"_test = [] {
        const auto matrix = scalar_general_matrix<T>();
        const std::vector<T> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        const std::vector<T> b = matrix * expected;
        const auto solver = stable_biconjugate_gradient{matrix};
        const std::vector<T> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };

    "block"_test = [] {
        const auto block_matrix = block_general_matrix<T>();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = block_matrix * expected;
        const auto solver = stable_biconjugate_gradient{block_matrix};
        const std::vector<std::array<T, 2>> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };

    constexpr size_t N = 1000;

    // "random_matrix"_test = [] {
    //     const auto matrix = random_general_matrix(N);
    //     const std::vector<T> expected = random_vector(N);
    //     const std::vector<T> b = matrix * expected;
    //     const auto solver = stable_biconjugate_gradient{matrix};

    //     const std::vector<T> x = solver.solve(b);
    //     std::cerr << "Iterations count: " << solver.iterations() << std::endl;

    //     // const T diff = norm<Inf>(x - expected);
    //     // expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    // };

    // "random_matrix_with_ilu0_preconditioner"_test = [] {
    //     const auto matrix = random_general_matrix(N);
    //     const std::vector<T> expected = random_vector(N);
    //     const std::vector<T> b = matrix * expected;
    //     auto solver = stable_biconjugate_gradient{matrix};

    //     solver.preconditioner(std::make_unique<ilu0_preconditioner<T, uint32_t, size_t>>(sparse_matrix<T>{matrix}));

    //     const std::vector<T> x = solver.solve(b);
    //     std::cerr << "Iterations count: " << solver.iterations() << std::endl;

    //     // const T diff = norm<Inf>(x - expected);
    //     // expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    // };

    // "random_block_matrix"_test = [] {
    //     const auto block_matrix = random_block_general_matrix(N);
    //     const std::vector<std::array<T, 2>> expected = random_block_vector(N);
    //     const std::vector<std::array<T, 2>> b = block_matrix * expected;
    //     const auto solver = stable_biconjugate_gradient{block_matrix};

    //     const std::vector<std::array<T, 2>> x = solver.solve(b);
    //     std::cerr << "Iterations count: " << solver.iterations() << std::endl;

    //     // const T diff = norm<Inf>(x - expected);
    //     // expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    // };

    // "random_block_matrix_with_ilu0_preconditioner"_test = [] {
    //     const auto block_matrix = random_block_general_matrix(N);
    //     const std::vector<std::array<T, 2>> expected = random_block_vector(N);
    //     const std::vector<std::array<T, 2>> b = block_matrix * expected;
    //     auto solver = stable_biconjugate_gradient{block_matrix};

    //     solver.preconditioner(std::make_unique<ilu0_preconditioner<square_matrix<T, 2>, uint32_t, size_t>>(sparse_matrix<square_matrix<T, 2>>{block_matrix}));

    //     const std::vector<std::array<T, 2>> x = solver.solve(b);
    //     std::cerr << "Iterations count: " << solver.iterations() << std::endl;

    //     // const T diff = norm<Inf>(x - expected);
    //     // expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    // };
};

}