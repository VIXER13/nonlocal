//#include <iostream>
//#include <petsc.h>
////#include "finite_element_solver_base.hpp"
//#include "heat_equation_solver.hpp"
//#include "influence_functions.hpp"
//
//void print_mesh(const mesh::mesh_2d<double>& msh) {
//    std::cout << "nodes:\n";
//    for(size_t n = 0; n < msh.nodes_count(); ++n)
//        std::cout << msh.node(n)[0] << ' ' << msh.node(n)[1] << '\n';
//    std::cout << std::endl;
//
//    std::cout << "elements:\n";
//    for(size_t el = 0; el < msh.elements_count(); ++el) {
//        const auto& e = msh.element_2d(el);
//        for(size_t i = 0; i < e->nodes_count(); ++i)
//            std::cout << msh.node_number(el, i) << ' ';
//        std::cout << '\n';
//    }
//    std::cout << std::endl;
//
//    std::cout << "boundary:\n";
//    for(size_t b = 0; b < msh.boundary_groups_count(); ++b) {
//        std::cout << "group" << b << ":\n";
//        for(size_t el = 0; el < msh.elements_count(b); ++el) {
//            const auto& e = msh.element_1d(b, el);
//            for(size_t i = 0; i < e->nodes_count(); ++i)
//                std::cout << msh.node_number(b, el, i) << ' ';
//            std::cout << '\n';
//        }
//    }
//}
//
////void print_solver_data(const nonlocal::finite_element_solver_base<double, int>& solver) {
////    std::cout << "quad shifts:\n";
////    for(size_t i = 0; i < solver._quad_shifts.size(); ++i)
////        std::cout << solver._quad_shifts[i] << ' ';
////    std::cout << std::endl;
////
////    std::cout << "quad coords:\n";
////    for(size_t i = 0; i < solver._quad_coords.size(); ++i)
////        std::cout << solver._quad_coords[i][0] << ' ' << solver._quad_coords[i][1] << '\n';
////    std::cout << std::endl;
////
////    std::cout << "jacobi matrices:\n";
////    for(size_t i = 0; i < solver._jacobi_matrices.size(); ++i)
////        std::cout << solver._jacobi_matrices[i][0] << ' ' << solver._jacobi_matrices[i][1] << ' '
////                  << solver._jacobi_matrices[i][2] << ' ' << solver._jacobi_matrices[i][3] << '\n';
////    std::cout << std::endl;
////
////    std::cout << "Neighbours:\n";
////    for(size_t i = 0; i < solver._elements_neighbors.size(); ++i) {
////        std::cout << i << " : ";
////        for(size_t j = 0; j < solver._elements_neighbors[i].size(); ++j)
////            std::cout << solver._elements_neighbors[i][j] << ' ';
////        std::cout << std::endl;
////    }
////}
//
//int main(int argc, char** argv) {
//    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
//
//    try {
//        PetscMPIInt size = -1, rank = -1;
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//        int maxthr = omp_get_max_threads() / size;
//        omp_set_num_threads(maxthr ? maxthr : 1);
//
//        mesh::mesh_2d<PetscScalar, PetscInt> msh{argv[1]};
//
//        static constexpr double r = 0.05, p1 = 0.5;
//        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
//        nonlocal::heat::heat_equation_solver<double, int> fem_sol{msh};
//
//        const auto T = fem_sol.stationary(
//            { // Граничные условия
//                {   // Down
//                        [](const std::array<double, 2>&) { return 0; },
//                        nonlocal::heat::boundary_t::FLOW
//                },
//
//                {   // Up
//                        [](const std::array<double, 2>&) { return 0; },
//                        nonlocal::heat::boundary_t::FLOW
//                },
//
//                {   // Left
//                        [](const std::array<double, 2>&) { return 0; },
//                        nonlocal::heat::boundary_t::FLOW
//                },
//
//                {   // Right
//                        [](const std::array<double, 2>&) { return 0; },
//                        nonlocal::heat::boundary_t::FLOW
//                }
//            },
//            [](const std::array<double, 2>&) { return 0; }, // Правая часть
//            r,  // Радиус влияния
//            p1, // Вес
//            bell // Функция влияния
//        );
//
//        //nonlocal::heat::heat_equation_solver<PetscScalar, PetscInt> solver{msh};
//        //solver.test();
//
////        //solver.find_neighbors(0.51);
////        for(int i = 0; i < size; ++i) {
////            if (i == rank) {
////                std::cout << "RANK: " << rank << std::endl;
////                //print_mesh(msh);
////                //print_solver_data(solver);
////                std::cout << std::endl << std::endl;
////            }
////            MPI_Barrier(MPI_COMM_WORLD);
////        }
//    } catch(const std::exception& e) {
//        std::cerr << e.what() << std::endl;
//        return EXIT_FAILURE;
//    } catch(...) {
//        std::cerr << "Unknown error." << std::endl;
//        return EXIT_FAILURE;
//    }
//
//    return PetscFinalize();
//}


//static char help[] = "Testing MatCreateMPIAIJSumSeqAIJ().\n\n";
//
//#include <petscmat.h>
//#include "../Eigen/Eigen/Sparse"
//#include <iostream>
//#include <vector>
//#include <cstdlib>
//
//int main(int argc,char **argv)
//{
//    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, help); if (ierr) return ierr;
//
//    PetscMPIInt size = 0, rank = 0;
//    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
//    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
//
//    const size_t glob_size = 16000000,
//                 loc_size  = glob_size / size;
//    Eigen::SparseMatrix<double, Eigen::RowMajor> C(glob_size, glob_size);
//
//    const size_t all_count = 160000;
//    const size_t count = all_count / size;
//    std::vector<Eigen::Triplet<double>> triplets(count);
//
//    for(int i = 0; i < size; ++i) {
//        if (i == rank)
//            std::cout << rank * loc_size << ' ' << (rank+1) * loc_size << std::endl;
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    for(size_t i = 0; i < count; ++i)
//        //triplets[i] = Eigen::Triplet<double>{i % loc_size + rank * loc_size, i % glob_size, 5.2};
//        //triplets[i] = Eigen::Triplet<double>{std::rand() % glob_size, std::rand() % glob_size, 5.2};
//        triplets[i] = Eigen::Triplet<double>{std::rand() % loc_size + rank * loc_size, std::rand() % glob_size, 5.2};
//
//    C.setFromTriplets(triplets.cbegin(), triplets.cend());
//
//    PetscLogDouble start_time = 0;
//    PetscTime(&start_time);
//
//    Mat A, B;
//    ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, C.rows(), C.cols(), C.outerIndexPtr(), C.innerIndexPtr(), C.valuePtr(), &A);CHKERRQ(ierr);
//    ierr = MatCreateMPIAIJSumSeqAIJ(PETSC_COMM_WORLD, A, PETSC_DECIDE, PETSC_DECIDE, MAT_INITIAL_MATRIX, &B);CHKERRQ(ierr);
//
//    PetscLogDouble end_time = 0;
//    PetscTime(&end_time);
//
//    if (rank == 0)
//        std::cout << end_time - start_time << std::endl;
//
//    PetscInt first = 0, last = 0;
//    ierr = MatGetOwnershipRange(B, &first, &last);CHKERRQ(ierr);
//
////    for(int i = 0; i < size; ++i) {
////        if (i == rank)
////            std::cout << first << ' ' << last << std::endl;
////        MPI_Barrier(MPI_COMM_WORLD);
////    }
//    //ierr = MatView(B, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//    ierr = MatDestroy(&B);CHKERRQ(ierr);
//    ierr = MatDestroy(&A);CHKERRQ(ierr);
//    return PetscFinalize();
//}







//static char help[] ="Tests MatPtAP() \n";
//
//#include <petscmat.h>
//#include <cstdlib>
//#include <iostream>
//#include "../Eigen/Eigen/Sparse"
//
//int main(int argc,char **argv)
//{
//    PetscErrorCode ierr;
//    Mat            A;
//    PetscInt       i1[] = {0, 3, 5}, i2[] = {0,0,0};//i2[] = {0,2,5};
//    PetscInt       j1[] = {0, 1, 3, 1, 2}, j2[] = {0, 0, 0, 0, 0};//j2[] = {0, 2, 1, 2, 3};
//    PetscScalar    a1[] = {1, 2, 4, 1, 2}, a2[] = {0, 0, 0, 0, 0};//a2[] = {2, 4, 1, 2, 1};
//    //PetscInt       pi1[] = {0,1,3}, pi2[] = {0,1,2};
//    //PetscInt       pj1[] = {0, 0, 1}, pj2[] = {1,0};
//    //PetscScalar    pa1[] = {1, 0.3, 0.5}, pa2[] = {0.8, 0.9};
//    MPI_Comm       comm;
//    PetscMPIInt    rank,size;
//
//
//
//
//    //C.outerIndexPtr(), C.innerIndexPtr()
//
//
//    ierr = PetscInitialize(&argc,&argv,NULL,help);if (ierr) return ierr;
//    comm = PETSC_COMM_WORLD;
//    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
//    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
//
//    int glob_size = 80000, loc_size = glob_size / size;
//    Eigen::SparseMatrix<double, Eigen::RowMajor> C{loc_size, glob_size};
//    int triplets_count = 16000000;
//    std::vector<Eigen::Triplet<double>> triplets(triplets_count / size);
//
//    for(size_t i = 0; i < triplets.size(); ++i)
//        triplets[i] = Eigen::Triplet<double>{std::rand() % loc_size, std::rand() % glob_size, 1.};
//
//    //if (size != 2)
//    //    SETERRQ(comm,PETSC_ERR_ARG_INCOMP,"You have to use two processor cores to run this example \n");
//    PetscLogDouble start_time = 0;
//    PetscTime(&start_time);
//    ierr = MatCreateMPIAIJWithArrays(comm, C.rows(), C.cols(), PETSC_DETERMINE, PETSC_DETERMINE, C.outerIndexPtr(), C.innerIndexPtr(), C.valuePtr(), &A);CHKERRQ(ierr);
//    PetscLogDouble end_time = 0;
//    PetscTime(&end_time);
//
//    std::cout << end_time - start_time << std::endl;
//    //ierr = MatView(A,NULL);CHKERRQ(ierr);
//    ierr = MatDestroy(&A);CHKERRQ(ierr);
//    ierr = PetscFinalize();
//    return ierr;
//}