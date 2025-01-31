#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>


/**
 * @brief Class representing an operator in a quantum system.
 * 
 * This class provides methods for initializing, manipulating, and diagonalizing operators represented as sparse matrices.
 */
class Operator {
private:

// INITIALIZATION :

    Eigen::SparseMatrix<double> O;
    int D; // dimension of the Hilbert space of the system
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
     * @brief Create the operator Identity.
     * 
     * @param size The size of the identity matrix.
     * @return Operator The identity operator.
     */
    static Operator Identity(int size);

    /**
     * @brief Add a matrix to an operand of type SparseMatrix with same size.
     * 
     * @param operand The matrix to add.
     * @return Operator The result of the addition.
     */
    Operator& operator + (const Operator& operand);

    /**
     * @brief Multiply a sparse matrix by a multiplicand of type SparseMatrix with same size.
     * 
     * @param multiplicand The matrix to multiply.
     * @return Operator The result of the multiplication.
     */
    Operator& operator * (const Operator& multiplicand);

    /**
     * @brief Multiply a sparse matrix by a vector with concomitant size.
     * 
     * @param vector The vector to multiply.
     * @return Eigen::Matrix<T, Eigen::Dynamic, 1> The result of the multiplication.
     */
    Eigen::VectorXd operator * (const Eigen::VectorXd& vector) const;

    /**
     * @brief Multiply a sparse matrix by a scalar.
     * 
     * @param scalar The scalar to multiply.
     * @return Operator The result of the multiplication.
     */
    Operator& operator * (double scalar);

// DIAGONALIZATION : 

    /**
     * @brief Calculate the approximate eigenvalues and eigenvectors of the Hamiltonian using the Implicitly Restarted Lanczos Method.
     * 
     * @param nb_eigen The number of eigenvalues to calculate.
     * @param eigenvectors An empty matrix to store the eigenvectors.
     * @return Eigen::Matrix<double> The vector of eigenvalues.
     */
    Eigen::VectorXcd IRLM_eigen(int nb_eigen) const;

    /**
     * @brief Calculate the approximate eigenvalues and eigenvectors of the Hamiltonian using the Full Orthogonalization Lanczos Method.
     * 
     * @param nb_iter The number of iterations.
     * @param v_0 The initial vector.
     * @param V An empty matrix to store the Ritz vectors.
     * @param eigenvectors An empty matrix to store the eigenvectors.
     * @return Eigen::Matrix<double> The vector of eigenvalues.
     */
    Eigen::VectorXd FOLM_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& eigenvectors) const;

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
    double gap_ratio();

    /** 
    * @brief Add a potential term to the operator.
    */
    void add_chemical_potential(double mu, int n);

    /**
    * @brief Add an interaction term to the operator.
    */
    void add_interaction(double U, const Eigen::VectorXd& basis);

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
    * @param dmu The difference of chemical potential.
    * @param n The number of bosons.
    * @return double The boson density.
    * @warning This function modifies the values of the operator from mu to mu + dmu.
    */
    double boson_density(double dmu, int n);

    /**
     * @brief Calculate the isothermal compressibility of the system.
     *
     * @param dmu The difference of chemical potential.
     * @param n The number of bosons.
     * @return double The compressibility.
     * @warning This function modifies the values of the operator from mu to mu + dmu.
     */
    double compressibility(double dmu, int n);

};