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

template<bool Is_Symmetric, int UpLo = Eigen::Upper, class T, int Major, class I>
auto find_eigen_value(const Eigen::SparseMatrix<T, Major, I>& matrix, const Spectra::SortRule rule) {
    if constexpr (Is_Symmetric) {
        Spectra::SparseSymMatProd<T, UpLo, Major, I> op{matrix};
        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<T, UpLo, Major, I>> eigs{op, 1, 256};
        eigs.init();
        eigs.compute(rule);
        return eigs.eigenvalues();
    } else {
        Spectra::SparseGenMatProd<T, Major, I> op{matrix};
        Spectra::GenEigsSolver<Spectra::SparseGenMatProd<T, Major, I>> eigs{op, 1, 256};
        eigs.init();
        eigs.compute(rule);
        return eigs.eigenvalues();
    }
}

}

#endif