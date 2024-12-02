#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <stdexcept>

class Hamiltonian {
private:
    Eigen::SparseMatrix<double> smat;
    double ref; // threshold under which a value is considered null
    long long m, n, D; // number of sites, bosons, and Hilbert space dimension

    long long factorial(long long n); // calculate factorial n with a recursive function
    long long dimension(long long m, long long n); // calculate the dimension of the Hilbert space for n bosons on m sites

    /* implement the Lanczos algorithm for a sparse matrix for nb_iter iterations starting with vector v_0 */
    void lanczos_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V);

    void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors); // sort eigenvalues and eigenvectors in descending order

public:
    Hamiltonian(long long m, long long n, double threshold = 1e-6); // constructor

    /* add a matrix to an operand of type SparseMatrix with same size */
    Hamiltonian operator + (const Hamiltonian& operand);

    /* multiply a sparse matrix by a multiplicand of type SparseMatrix with same size */
    Hamiltonian operator * (const Hamiltonian& multiplicand);

    /* multiply a sparse matrix by a vector with concomitant size */
    Eigen::VectorXd operator * (const Eigen::VectorXd& vector);

    /* calculate the approximate eigenvalues and eigenvectors of the Hamiltonian in Krylov space using the Lanczos algorithm */
    Eigen::VectorXd lanczos_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V);

    /* calculate the exact eigenvalues and eigenvectors of the Hamiltonian by an exact diagonalization */
    Eigen::VectorXd exact_eigen(Eigen::MatrixXd& eigenvectors);
};

#endif // HAMILTONIAN_H
