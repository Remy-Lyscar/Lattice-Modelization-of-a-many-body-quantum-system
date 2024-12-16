#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "operator.h"


// OPERATOR CONSTRUCTOR :
Operator::Operator(Eigen::SparseMatrix<double>&& smatrix) : smat(std::move(smatrix)), D(smatrix.rows()), ref(1e-6) {}


//UTILITY FUNCTIONS :
int Operator::size() const {
	return D;
}

// BASIS OPERATIONS : 

/* add a matrix to an operand of type SparseMatrix with same size */
Operator Operator::operator + (const Operator& operand) const {
    if (this->smat.rows() != operand.smat.rows() || this->smat.cols() != operand.smat.cols()) { // verify that the operands have matching size
        throw std::invalid_argument("Matrix should have matching size.");
    }
    Eigen::SparseMatrix<double> smatrix(this->smat.rows(), this->smat.cols());
    Operator result(std::move(smatrix));
    result.smat = (this->smat + operand.smat).pruned(ref); // removes elements smaller than ref
    return result;
}

/* multiply a sparse matrix by a multiplicand of type SparseMatrix with same size */
Operator Operator::operator * (const Operator& multiplicand) const {
    if (this->smat.cols() != multiplicand.smat.rows()) {
        throw std::invalid_argument("Number of columns of multiplier must equal number of rows of multiplicand."); // verify that the number of columns of multiplier equals number of rows of multiplicand
    }
    Eigen::SparseMatrix<double> smatrix(this->smat.rows(), multiplicand.smat.cols());
    Operator result(std::move(smatrix));
    result.smat = (this->smat * multiplicand.smat).pruned(ref); // removes elements smaller than ref
    return result;
}

/* multiply a sparse matrix by a vector with concomitant size */
Eigen::VectorXd Operator::operator * (const Eigen::VectorXd& vector) const {
    if (this->smat.cols() != vector.size()) { // verify that number of matrix columns equals vector size
        throw std::invalid_argument("Number of matrix columns must equal vector size.");
    }
    return this->smat * vector;
}

// DIAGONALIZATION : 

/* implement the Lanczos algorithm for a sparse matrix for nb_iter iterations starting with vector v_0 */
void Operator::lanczos_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V) const {
    T.resize(nb_iter, nb_iter); //resize the matrix T to a matching size
    V.resize(D, nb_iter); //resize the matrix V to a matching size

    v_0.normalize(); // normalize the starting vector v_0
    V.col(0) = v_0; // put the first vector v_0 in the matrix V

    double alpha, beta;

    //Lanczos algorithm for nb_iter iterations
    for (int i = 0; i < nb_iter; i++) {
        Eigen::VectorXd w = this->smat * V.col(i); // calculate the next vector w 
        alpha = (V.col(i)).dot(w);
        for (int j = 0; j < i; j++) {
            w = w - (V.col(j)).dot(w) * V.col(j); // orthogonalize the vector w with respect to the previous vectors of the Krylov basis
        }
        beta = w.norm(); // calculate the norm of the vector w
        if (beta < ref) {
            break; // if beta is null or almost null the algorithm stops
        }
        else {
            w = w / beta; // normalize the vector w
            if (i + 1 < nb_iter) {
                V.col(i + 1) = w; // add the vector w to the matrix V of vectors of the Krylov basis
            }
            T(i, i) = alpha; // add the ith diagonal element of the tridiagonal matrix T
            if (i > 0) {
                T(i, i - 1) = beta; // add the ith non-diagonal element of the tridiagonal matrix T
                T(i - 1, i) = beta; // add the ith non-diagonal element of the tridiagonal matrix T
            }
        }
    }
}

/* sort the eigenvalues and the eigenvectors by the eigenvalue in descending order */
void Operator::sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) const {
    if (eigenvalues.size() != eigenvectors.cols()) {
        throw std::invalid_argument("Eigenvalues and eigenvectors must have matching sizes.");
    }
    for (int i = 0; i < eigenvalues.size(); i++) {
        for (int j = i + 1; j < eigenvalues.size(); j++) {
            if (eigenvalues[i] < eigenvalues[j]) { // change i and j for ascending order
                std::swap(eigenvalues[i], eigenvalues[j]);
                eigenvectors.col(i).swap(eigenvectors.col(j));
            }
        }
    }
}

/* calculate the approximate eigenvalues and eigenvectors of the hamiltonian in Krylov space using the Lanczos algorithm */
Eigen::VectorXd Operator::lanczos_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V, Eigen::MatrixXd& eigenvectors) const {
    Eigen::MatrixXd T(nb_iter,nb_iter); // initialize the tridiagonal matrix T 
    lanczos_diag(nb_iter, v_0, T, V); // tridiagonalize the hamiltonian using Lanczos algorithm

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(T); // solve the eigen problem for T
    if (eigensolver.info() != Eigen::Success) { // verify if the eigen search is a success
        throw std::runtime_error("Eigenvalue computation failed.");
    }
    eigenvectors = V * eigensolver.eigenvectors(); // eigenvectors of the hamiltonian
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues(); // eigenvalues of the hamiltonian
    sort_eigen(eigenvalues, eigenvectors);
    return eigenvalues;
}

/* calculate the exact eigenvalues and eigenvectors of the hamiltonian by an exact diagonalization */
Eigen::VectorXd Operator::exact_eigen(Eigen::MatrixXd& eigenvectors) const {
    Eigen::MatrixXd dense_smat = Eigen::MatrixXd(this->smat); // convert sparse matrix to dense matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(dense_smat); // solve the eigen problem for the hamiltonian
    if (eigensolver.info() != Eigen::Success) { // verify if the eigen search is a success
        throw std::runtime_error("Eigenvalue computation failed.");
    }
    eigenvectors = eigensolver.eigenvectors(); // eigenvectors of the hamiltonian
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues(); // eigenvalues of the hamiltonian
    sort_eigen(eigenvalues, eigenvectors);
    return eigenvalues;
}

// THERMODYNAMICAL FUNCTIONS :

/* calculate the partition function Z for an ALREADY diagonalized hamiltonian */
double Operator::partition_function(const Eigen::VectorXd& eigenvalues, double temperature) const {
    const double k_B = 1.380649e-23; // Boltzmann constant
    double beta = 1 / (k_B * temperature);
    double Z = 0;
	for (int i = 0; i < eigenvalues.size(); i++) { // calculate the partition function using properties of the trace
        Z += exp(-beta * eigenvalues[i]); 
    }
    return Z;
}

/* calculate the canonical density matrix for an ALREADY diagonalized hamiltonian */
void Operator::canonical_density_matrix(const Eigen::VectorXd& eigenvalues, double temperature) const {
    throw std::logic_error("This function has not been implemented yet.");
}
