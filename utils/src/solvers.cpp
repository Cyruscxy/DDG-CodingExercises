#include "solvers.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Compute the inverse of a sparse diagonal matrix.
 *
 * Input: A sparse diagonal matrix <M>.
 * Returns: The inverse of M, which is also a sparse diagonal matrix.
 */
SparseMatrix<double> sparseInverseDiagonal(SparseMatrix<double>& M) {

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList(M.rows());
    SparseMatrix<double> inv(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeffRef(i, i)));
    }
    inv.setFromTriplets(tripletList.begin(), tripletList.end());
    return inv;
}

/*
 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
 *
 * Input: <A>, the complex sparse matrix whose eigendecomposition is being computed; and <x>, the current guess for the
 * smallest eigenvector
 * Returns: The residual
 */
double residual(const SparseMatrix<std::complex<double>>& A, const Vector<std::complex<double>>& x) {

    Vector<std::complex<double>> x_ht(x.size());
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
        x_ht[i] = conj(x[i]);
    }

    const auto lambda = (x_ht.transpose() * A * x).sum() / (x_ht.transpose() * x).sum();

    return (A * x - lambda * x).norm() / x.norm();
}

/*
 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A, and x is the corresponding eigenvector.
 *
 * Input: <A>, the complex positive definite sparse matrix whose eigendecomposition is being computed.
 * Returns: The smallest eigenvector of A.
 */
Vector<std::complex<double>> solveInversePowerMethod(const SparseMatrix<std::complex<double>>& A) {
    auto L = A;
    Vector<std::complex<double>> y(A.cols());
    y.setRandom();

    while (residual(A, y) > 1e-10) {
		y -= Vector<std::complex<double>>::Constant(y.rows(), y.sum() / std::complex<double>(y.rows(), 0)) ;
        y = geometrycentral::solvePositiveDefinite<std::complex<double>>(L, y);
        y /= y.norm();
    }

    return y;
}