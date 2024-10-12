#ifndef NONLOCAL_EIGEN_VALUES_HPP
#define NONLOCAL_EIGEN_VALUES_HPP

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <iostream>
#include <omp.h>

namespace nonlocal::slae {

template<int UpLo = Eigen::Upper, class T, int Major, class I>
T find_eigen_value(const Eigen::SparseMatrix<T, Major, I>& matrix, const bool is_symmetric, const Spectra::SortRule rule) {
    if (is_symmetric) {
        Spectra::SparseSymMatProd<T, UpLo, Major, I> op{matrix};
        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<T, UpLo, Major, I>> eigs{op, 1, 256};
        eigs.init();
        eigs.compute(rule);
        return eigs[1];
    } else {
        Spectra::SparseGenMatProd<T, Major, I> op{mmm};
        Spectra::GenEigsSolver<Spectra::SparseGenMatProd<T, Major, I>> eigs{op, 1, 256};
        eigs.init();
        eigs.compute(rule);
        return eigs[1];
    }
}

}

#endif