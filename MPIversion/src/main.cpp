#include <iostream>
#include <petsc.h>
#include "finite_element_solver_base.hpp"

void print_mesh(const mesh::mesh_2d<double>& msh) {
    std::cout << "nodes:\n";
    for(size_t n = 0; n < msh.nodes_count(); ++n)
        std::cout << msh.node(n)[0] << ' ' << msh.node(n)[1] << '\n';
    std::cout << std::endl;

    std::cout << "elements:\n";
    for(size_t el = 0; el < msh.elements_count(); ++el) {
        const auto& e = msh.element_2d(el);
        for(size_t i = 0; i < e->nodes_count(); ++i)
            std::cout << msh.node_number(el, i) << ' ';
        std::cout << '\n';
    }
    std::cout << std::endl;

    std::cout << "boundary:\n";
    for(size_t b = 0; b < msh.boundary_groups_count(); ++b) {
        std::cout << "group" << b << ":\n";
        for(size_t el = 0; el < msh.elements_count(b); ++el) {
            const auto& e = msh.element_1d(b, el);
            for(size_t i = 0; i < e->nodes_count(); ++i)
                std::cout << msh.node_number(b, el, i) << ' ';
            std::cout << '\n';
        }
    }
}

void print_solver_data(const nonlocal::finite_element_solver_base<double, int>& solver) {
    std::cout << "quad shifts:\n";
    for(size_t i = 0; i < solver._quad_shifts.size(); ++i)
        std::cout << solver._quad_shifts[i] << ' ';
    std::cout << std::endl;

    std::cout << "quad coords:\n";
    for(size_t i = 0; i < solver._quad_coords.size(); ++i)
        std::cout << solver._quad_coords[i][0] << ' ' << solver._quad_coords[i][1] << '\n';
    std::cout << std::endl;

    std::cout << "jacobi matrices:\n";
    for(size_t i = 0; i < solver._jacobi_matrices.size(); ++i)
        std::cout << solver._jacobi_matrices[i][0] << ' ' << solver._jacobi_matrices[i][1] << ' '
                  << solver._jacobi_matrices[i][2] << ' ' << solver._jacobi_matrices[i][3] << '\n';
    std::cout << std::endl;

    std::cout << "Neighbours:\n";
    for(size_t i = 0; i < solver._elements_neighbors.size(); ++i) {
        std::cout << i << " : ";
        for(size_t j = 0; j < solver._elements_neighbors[i].size(); ++j)
            std::cout << solver._elements_neighbors[i][j] << ' ';
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    try {
        PetscMPIInt size = -1, rank = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        mesh::mesh_2d<PetscScalar, PetscInt> msh{argv[1]};
        nonlocal::finite_element_solver_base<PetscScalar, PetscInt> solver{msh};
        solver.find_neighbors(0.51);
        for(int i = 0; i < size; ++i) {
            if (i == rank) {
                std::cout << "RANK: " << rank << std::endl;
                //print_mesh(msh);
                print_solver_data(solver);
                std::cout << std::endl << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
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