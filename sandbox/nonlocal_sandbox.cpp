#include <omp.h>

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include <solvers/slae/independent_symmetric_matrix_vector_product.hpp>

constexpr size_t ITERS = 100;

namespace Native_NS {
template <class T, class I>
void dot(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
         const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
         Eigen::Matrix<T, Eigen::Dynamic, 1>& result) {
    result.resize(matrix.rows());
    for (size_t i = 0; i < ITERS; ++i) {
        result.setZero();
        result = matrix.template selfadjointView<Eigen::Upper>() * vector;  // Встроенная оптимизированная операция
    }
}
}  // namespace Native_NS

namespace Iterator_NS {
template <class T, class I>
void dot(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
         const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
         Eigen::Matrix<T, Eigen::Dynamic, 1>& result) {
    result.resize(matrix.rows());
    const I n = static_cast<I>(matrix.rows());
    const int max_threads = omp_get_max_threads();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> local_results(
        static_cast<size_t>(max_threads),
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n));

    for (size_t i = 0; i < ITERS; ++i) {
        result.setZero();
#pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            auto& local = local_results[static_cast<size_t>(tid)];

#pragma omp for schedule(dynamic, 64)
            for (I i = 0; i < n; ++i) {
                for (typename Eigen::SparseMatrix<T, Eigen::RowMajor, I>::InnerIterator it(
                        matrix, i);
                    it; ++it) {
                    const I j = static_cast<I>(it.col());
                    const T a = it.value();
                    local(i) += a * vector(j);
                    if (j != i) {
                        local(j) += a * vector(i);
                    }
                }
            }
        }

        result = local_results.front();
        for (size_t t = 1; t < local_results.size(); ++t) {
            result.noalias() += local_results[t];
        }
    }
}
}  // namespace Iterator_NS

namespace NonLocal_NS {

template <class T, class I>
void dot(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
         const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
         Eigen::Matrix<T, Eigen::Dynamic, 1>& result) {
            result.resize(matrix.cols());
            nonlocal::slae::independent_symmetric_matrix_vector_product<T, I> product(matrix);
            for (size_t i = 0; i < ITERS; ++i) {
                product.matrix_vector_product(result, vector);
            }
         }
} // namespace NonLocal_NS

// Замер времени выполнения функции
template <typename Func>
double measure_time(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    return elapsed.count();
}

// Создание разреженной матрицы с заданной плотностью
template <class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I> create_sparse_matrix(
    I rows, I cols, double density = 0.01, unsigned int seed = 42) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> dist(0.0, 1.0);
    std::uniform_int_distribution<I> col_dist(0, cols - 1);

    Eigen::SparseMatrix<T, Eigen::RowMajor, I> matrix(rows, cols);
    std::vector<Eigen::Triplet<T, I>> tripletList;
    tripletList.reserve(rows * cols * density);

    for (I i = 0; i < rows; ++i) {
        I nnz_in_row = static_cast<I>(cols * density);
        for (I j = 0; j < nnz_in_row; ++j) {
            I col = col_dist(gen);
            T val = dist(gen);
            if (col >= i)
                tripletList.emplace_back(i, col, val);
        }
        tripletList.emplace_back(i, i, 1.0);
    }

    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    matrix.makeCompressed();
    return matrix;
}

// Создание случайного вектора
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> create_random_vector(
    Eigen::Index size, unsigned int seed = 43) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> dist(0.0, 1.0);
    Eigen::Matrix<T, Eigen::Dynamic, 1> vector(size);
    for (Eigen::Index i = 0; i < size; ++i) {
        vector(i) = dist(gen);
    }
    return vector;
}

// Проверка близости двух векторов
template <class T>
bool verify_results(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v1,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& v2,
                    T tolerance = 1e-6) {
    if (v1.size() != v2.size()) 
        return false;
    const T max = (v1 - v2).cwiseAbs().maxCoeff();
    std::cout << max << std::endl;
    return max < tolerance;
}

// ============================================================================
// Основная функция
// ============================================================================
int main(int argc, char* argv[]) {
    // Параметры
    int size = 10000;      // Размер матрицы и вектора
    double density = 0.1;  // Плотность ненулевых элементов (1%)
    int num_threads = 4;    // Количество потоков OpenMP

    // Чтение параметров из командной строки
    if (argc > 1) size = std::atoi(argv[1]);
    if (argc > 2) density = std::atof(argv[2]);
    if (argc > 3) num_threads = std::atoi(argv[3]);

    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);

    std::cout << "==============================================\n";
    std::cout << "Sparse Matrix-Vector Multiplication Benchmark\n";
    std::cout << "==============================================\n";
    std::cout << "Matrix size: " << size << " x " << size << "\n";
    std::cout << "Density: " << (density * 100) << "%\n";
    std::cout << "Threads: " << num_threads << "\n";
    std::cout << "==============================================\n\n";

    // Создание данных
    std::cout << "Creating sparse matrix... ";
    auto matrix = create_sparse_matrix<double, int>(size, size, density);
    std::cout << "Done (nnz: " << matrix.nonZeros() << ")\n";

    std::cout << "Creating random vector... ";
    auto vector = create_random_vector<double>(size);
    std::cout << "Done\n\n";

    // Результаты
    Eigen::VectorXd result_native, result_iterator, result_nonlocal;
    double time_native = 0.0, time_iterator = 0.0, time_nonlocal = 0.0;

    // Native_NS
    time_native = measure_time([&]() { Native_NS::dot(matrix, vector, result_native); });

    // Iterator_NS
    time_iterator = measure_time([&]() { Iterator_NS::dot(matrix, vector, result_iterator); });

    // NonLocal_NS
    time_nonlocal = measure_time([&]() { NonLocal_NS::dot(matrix, vector, result_nonlocal); });

    // Вывод результатов
    std::cout << "\n==============================================\n";
    std::cout << "Results:\n";
    std::cout << "==============================================\n";
    std::cout << "Native_NS time:    " << time_native << " ms\n";
    std::cout << "Iterator_NS time:  " << time_iterator << " ms\n";
    std::cout << "NonLocal_NS time:  " << time_nonlocal << " ms\n";

    std::cout << "Speedup (NonLocal/Iterator): " << time_nonlocal / time_iterator << "x\n";
    std::cout << "Speedup (NonLocal/Native): " << time_nonlocal / time_native << "x\n";

    // Проверка корректности
    std::cout << "\nVerification: ";
    if (verify_results(result_native, result_iterator) && verify_results(result_native, result_nonlocal)) {
        std::cout << "PASSED ✓\n";
    } else {
        std::cout << "FAILED ✗\n";
    }

    std::cout << "\n==============================================\n";

    return 0;
}