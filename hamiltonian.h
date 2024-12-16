#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <vector>

class Hamiltonian {
private:
    Eigen::SparseMatrix<double> mat;
    std::vector<std::vector<int>> neighbours;
    int m, n;
    double J, U, mu;

    // DIMENSION OF THE HILBERT SPACE

    /* calculate factorial n with a recursive function */
    int factorial(int n) const;

    /* calculate the dimension of the Hilbert space for n bosons on m sites */
    int dimension(int m, int n) const;

    // ELEMENTARY FUNCTIONS

    /* calculate the sum of the elements of a vector between 2 index */
    int sum(const Eigen::VectorXd& state, int index1, int index2) const;

    // INITIALIZE THE HILBERT SPACE BASIS

    /* calculate the next state of the Hilbert space in lexicographic order */
    bool next_lexicographic(Eigen::VectorXd& state, int m, int n) const;

    /* creates a matrix that has the vectors of the Hilbert space basis in columns */
    Eigen::MatrixXd init_lexicographic(int m, int n) const;

    // SORT THE HILBERT SPACE BASIS TO FACILITATE CALCULUS

    /* calculate the unique tag of the kth column of the matrix */
    double calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k) const;

    /* calculate and store the tags of each state of the Hilbert space basis */
    Eigen::VectorXd calculate_tags(const Eigen::MatrixXd& basis, const std::vector<int>& primes) const;

    /* sort the states of the Hilbert space by ascending order compared by their tags */
    void sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) const;

    // FILL THE HAMILTONIAN OF THE SYSTEM

    /* gives the index of the wanted tag x by the Newton method */
    int search_tag(const Eigen::VectorXd& tags, double x) const;

    /* fill the hopping term of the Hamiltonian */
    void fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, const std::vector<int>& primes, Eigen::SparseMatrix<double>& hmatrix, double J) const;

    /* fill the interaction term of the Hamiltonian */
    void fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double U) const;

    /* fill the chemical potential term of the Hamiltonian */
    void fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double mu) const;

public:
    /* constructor for the Hamiltonian class */
    Hamiltonian(const std::vector<std::vector<int>>& neighbours, int m, int n, double J, double U, double mu);

    /* get the Hamiltonian matrix */
    Eigen::SparseMatrix<double> getHamiltonian() const;
};

#endif // HAMILTONIAN_H
