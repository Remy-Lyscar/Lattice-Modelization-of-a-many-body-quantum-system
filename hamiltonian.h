#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include <cmath>

class Hamiltonian {
private:
    Eigen::SparseMatrix<double> mat; 
    std::vector<std::vector<int>> neighbours;
    int m, n; 
    double J, U, mu; 

    int factorial(int n);
    int dimension(int m, int n);
    int sum(const Eigen::VectorXd& state, int index1, int index2);
    bool next_lexicographic(Eigen::VectorXd& state, int m, int n);
    Eigen::MatrixXd init_lexicographic(int m, int n);
    double calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k);
    Eigen::VectorXd calculate_tags(const Eigen::MatrixXd& basis, const std::vector<int>& primes);
    void sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis);
    int search_tag(const Eigen::VectorXd& tags, double x);
    void fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, Eigen::SparseMatrix<double>& hmatrix, double J);
    void fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double U);
    void fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double mu);

public:
    Hamiltonian(const std::vector<std::vector<int>>& neighbours, int m, int n, double J, double U, double mu);
    Eigen::SparseMatrix<double> getHamiltonian();
};