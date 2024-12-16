#ifndef OPERATOR_H
#define OPERATOR_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <stdexcept>

class Operator {
private:
	// INITIALIZATION :
    Eigen::SparseMatrix<double> smat;
    int D;
    double ref; // threshold under which a value is considered null

    // DIAGONALIZATION : 

    /* implement the Lanczos algorithm for a sparse matrix for nb_iter iterations starting with vector v_0 */
    void lanczos_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V) const;

    /* sort eigenvalues and eigenvectors in descending order */
    void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) const;

public:
    // OPERATOR CONSTRUCTOR 
    Operator(Eigen::SparseMatrix<double>&& smatrix);

    //UTILITY FUNCTIONS :
    int size() const;

    // BASIS OPERATIONS : 

    /* add a matrix to an operand of type SparseMatrix with same size */
    Operator operator + (const Operator& operand) const;

    /* multiply a sparse matrix by a multiplicand of type SparseMatrix with same size */
    Operator operator * (const Operator& multiplicand) const;

    /* multiply a sparse matrix by a vector with concomitant size */
    Eigen::VectorXd operator * (const Eigen::VectorXd& vector) const;

    // DIAGONALIZATION : 

    /* calculate the approximate eigenvalues and eigenvectors of the Hamiltonian in Krylov space using the Lanczos algorithm */
    Eigen::VectorXd lanczos_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V, Eigen::MatrixXd& eigenvectors) const;

    /* calculate the exact eigenvalues and eigenvectors of the Hamiltonian by an exact diagonalization */
    Eigen::VectorXd exact_eigen(Eigen::MatrixXd& eigenvectors) const;

    // THERMODYNAMICAL FUNCTIONS :

    /* calculate the partition function Z for an ALREADY diagonalized Hamiltonian */
    double partition_function(const Eigen::VectorXd& eigenvalues, double temperature) const;

    void canonical_density_matrix(const Eigen::VectorXd& eigenvalues, double temperature) const;
};

#endif // OPERATOR_H
