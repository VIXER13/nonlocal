#include <petsc.h>
#include "influence_functions.hpp"
#include "heat_equation_solver.hpp"

namespace {

void save_raw_data(const std::shared_ptr<mesh::mesh_2d<double>>& mesh, const Eigen::Matrix<double, Eigen::Dynamic, 1>& T) {
    std::ofstream Tout{"T.csv"};
    for(size_t i = 0; i < mesh->nodes_count(); ++i)
        Tout << mesh->node(i)[0] << ',' << mesh->node(i)[1] << ',' << T[i] << '\n';
}

}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);
        omp_set_num_threads(1);

        static constexpr double r = 0.2, p1 = 1;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);

        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_info = std::make_shared<mesh::mesh_info<double, int>>(mesh);
        mesh_info->find_neighbours(r);
        nonlocal::heat::heat_equation_solver<double, int> fem_sol{mesh_info};

        const auto T = fem_sol.stationary(
            { // Граничные условия
                {   // Down
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                },

                {   // Right
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                },

                {   // Up
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                },

                {   // Left
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return x[0]*x[0] + x[1]*x[1]; },
                }
            },
            [](const std::array<double, 2>&) { return -4; }, // Правая часть
            p1, // Вес
            bell // Функция влияния
        );
//
        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            fem_sol.save_as_vtk("heat.vtk", T);
            save_raw_data(mesh, T);
            std::cout << fem_sol.integrate_solution(T) << std::endl;
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return PetscFinalize();
}



//static char help[] = "Tests MatCreateMPISBAIJWithArrays().\n\n";
//
//#include <petscmat.h>
//#include "../../Eigen/Eigen/Sparse"
//
//int main(int argc,char **args)
//{
//    Mat            A;
//    PetscErrorCode ierr;
//    PetscInt       m = 4, bs = 1,ii[5],jj[7];
//    PetscMPIInt    size,rank;
//    PetscScalar    aa[7];
//
//    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
//    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
//    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//    if (size != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only for two processes");
//
//    Eigen::SparseMatrix<double, Eigen::RowMajor> C{4, 8};
//
//    if (!rank) {
//        ii[0] = 0; ii[1] = 2; ii[2] = 5; ii[3] = 7; ii[4] = 7;
//        jj[0] = 0; jj[1] = 1; jj[2] = 1; jj[3] = 2; jj[4] = 6; jj[5] = 3; jj[6] = 7;
//        aa[0] = 0; aa[1] = 1; aa[2] = 2; aa[3] = 3; aa[4] = 4; aa[5] = 5; aa[6] = 6;
//        /*  0 1
//              1  2       6
//                 3          7 */
//
//        C.coeffRef(0, 0) = 0;
//        C.coeffRef(0, 1) = 1;
//        C.coeffRef(1, 1) = 2;
//        C.coeffRef(1, 2) = 3;
//        C.coeffRef(1, 6) = 4;
//        C.coeffRef(2, 3) = 5;
//        C.coeffRef(2, 7) = 5;
//    } else {
//        ii[0] = 0; ii[1] = 2; ii[2] = 4; ii[3] = 6; ii[4] = 7;
//        jj[0] = 4; jj[1] = 5; jj[2] = 5; jj[3] = 7; jj[4] = 6; jj[5] = 7; jj[6] = 7;
//        aa[0] = 8; aa[1] = 9; aa[2] = 10; aa[3] = 11; aa[4] = 12; aa[5] = 13; aa[6] = 14;
//        /*    4  5
//                 5   7
//                   6 7
//                     7 */
//
//        C.coeffRef(0, 4) = 8;
//        C.coeffRef(0, 5) = 9;
//        C.coeffRef(1, 5) = 10;
//        C.coeffRef(1, 7) = 11;
//        C.coeffRef(2, 6) = 12;
//        C.coeffRef(2, 7) = 13;
//        C.coeffRef(3, 7) = 14;
//    }
//    C.makeCompressed();
//
//    ierr = MatCreateMPISBAIJWithArrays(PETSC_COMM_WORLD, bs, C.rows(), PETSC_DECIDE, PETSC_DECIDE, C.cols(),
//                                       C.outerIndexPtr(), C.innerIndexPtr(), C.valuePtr(), &A);CHKERRQ(ierr);
//    //ierr = MatCreateMPISBAIJWithArrays(PETSC_COMM_WORLD,bs,m,m,PETSC_DECIDE,PETSC_DECIDE,ii,jj,aa,&A);CHKERRQ(ierr);
//    //ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,m,m,PETSC_DECIDE,PETSC_DECIDE,ii,jj,aa,&A);CHKERRQ(ierr);
//    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//    ierr = MatDestroy(&A);CHKERRQ(ierr);
//    ierr = PetscFinalize();
//    return ierr;
//}