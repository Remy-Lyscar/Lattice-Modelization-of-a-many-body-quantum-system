#ifndef OPERATOR_H
#define OPERATOR_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>



/**
 * @brief Class representing an operator in a quantum system.
 * 
 * This class provides methods for initializing, manipulating, and diagonalizing operators represented as sparse matrices.
 */
class Operator {
private:

// INITIALIZATION :

    Eigen::SparseMatrix<double> O;
    int D;
    double ref; // threshold under which a value is considered null

// DIAGONALIZATION : 

    /* implement the Full Orthogonalization Lanczos Method for a sparse matrix for nb_iter iterations starting with vector v_0 */
    void FOLM_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V) const;

    /* sort eigenvalues and eigenvectors in descending order */
    void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) const;

public:

// CONSTRUCTOR :

    /**
     * @brief Constructor for the Operator class.
     * 
     * @param smatrix The sparse matrix to initialize the operator.
     */
    Operator(Eigen::SparseMatrix<double>&& smatrix);

//UTILITY FUNCTIONS :

    /**
     * @brief Get the size of the matrix.
     * 
     * @return int The size of the matrix.
     */
    int size() const;

// BASIS OPERATIONS : 

    /**
     * @brief Add a matrix to an operand of type SparseMatrix with same size.
     * 
     * @param operand The matrix to add.
     * @return Operator The result of the addition.
     */
    Operator operator + (const Operator& operand) const;

    /**
     * @brief Multiply a sparse matrix by a multiplicand of type SparseMatrix with same size.
     * 
     * @param multiplicand The matrix to multiply.
     * @return Operator The result of the multiplication.
     */
    Operator operator * (const Operator& multiplicand) const;

    /**
     * @brief Multiply a sparse matrix by a vector with concomitant size.
     * 
     * @param vector The vector to multiply.
     * @return Eigen::Matrix<T, Eigen::Dynamic, 1> The result of the multiplication.
     */
    Eigen::VectorXd operator * (const Eigen::VectorXd& vector) const;

// DIAGONALIZATION : 

    /**
     * @brief Calculate the approximate eigenvalues and eigenvectors of the Hamiltonian using the Implicitly Restarted Lanczos Method.
     * 
     * @param k The number of eigenvalues to calculate.
     * @param eigenvectors An empty matrix to store the eigenvectors.
     * @return Eigen::Matrix<double> The vector of eigenvalues.
     */
    Eigen::VectorXd IRLM_eigen(int k, Eigen::MatrixXd& eigenvectors) const;

    /**
     * @brief Calculate the approximate eigenvalues and eigenvectors of the Hamiltonian using the Full Orthogonalization Lanczos Method.
     * 
     * @param nb_iter The number of iterations.
     * @param v_0 The initial vector.
     * @param V An empty matrix to store the Ritz vectors.
     * @param eigenvectors An empty matrix to store the eigenvectors.
     * @return Eigen::Matrix<double> The vector of eigenvalues.
     */
    Eigen::VectorXd FOLM_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V, Eigen::MatrixXd& eigenvectors) const;

    /**
     * @brief Calculate the exact eigenvalues and eigenvectors of the Hamiltonian by an exact diagonalization.
     * 
     * @param eigenvectors An empty matrix to store the eigenvectors.
     * @return Eigen::Matrix<double> The vector of eigenvalues.
     */
    Eigen::VectorXd exact_eigen(Eigen::MatrixXd& eigenvectors) const;


// PHASE TRANSITION CALCULATIONS :

    /**
     * @brief Calculate the order parameter of the system.
     * 
     * @param eigenvalues The vector of eigenvalues.
     * @param eigenvectors The matrix of eigenvectors.
     * @return double The order parameter.
     */
    double order_parameter(const Eigen::VectorXd& eigenvalues, const Eigen::MatrixXd& eigenvectors) const;

    /**
    * @brief Calculate the energy gap ratio of the system.
    *
    * @param eigenvalues The vector of eigenvalues.
    * @return double The energy gap ratio.
    */
    double gap_ratio(const Eigen::VectorXd& eigenvalues) const;

// THERMODYNAMICAL FUNCTIONS :

    /**
     * @brief Calculate the partition function Z for an already diagonalized Hamiltonian.
     * 
     * @param eigenvalues The vector of eigenvalues.
     * @param temperature The temperature.
     * @return double The partition function.
     */
    double partition_function(const Eigen::VectorXd& eigenvalues, double temperature) const;

    /**
     * @brief Calculate the canonical density matrix for an already diagonalized Hamiltonian.
     * 
     * @param eigenvalues The vector of eigenvalues.
     * @param temperature The temperature.
     */
    void canonical_density_matrix(const Eigen::VectorXd& eigenvalues, double temperature) const;

    /**
     * @brief Calculate the boson density of the system.
     * 
     * @param eigenvalues1 The vector of eigenvalues of the first Hamiltonian with mu.
     * @param eigenvalues2 The vector of eigenvalues of the second Hamiltonian with mu + dmu.
     * @param dmu The difference of chemical potential.
     * @return double The boson density.
     */
    double boson_density(const Eigen::VectorXd& eigenvalues1, const Eigen::VectorXd& eigenvalues2, double dmu ) const;

    /**
     * @brief Calculate the double compressibility of the system.
     * 
     * @param density1 The density of the first Hamiltonian with mu.
     * @param density2 The density of the second Hamiltonian with mu + dmu.
     * @param dmu The difference of chemical potential.
     * @return double The double compressibility.
     */
    double compressibility(double density1, double density2, double dmu) const;

};

#endif // OPERATOR_H
