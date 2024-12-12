#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "hamiltonian.h"
#include "operator.h"

int main() {
	int m = 3;
	int n = 3;

	double J = 1;
	double U = 1;
	double mu = 1;

	std::vector < std::vector<int>> neighbours = { {1}, {0,2}, {1} };

	Hamiltonian hmatrix (neighbours, m, n, J, U, mu);

	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();

	Operator H(std::move(smatrix));
	Eigen::MatrixXd eigenvectors;
	Eigen::VectorXd eigenvalues = H.exact_eigen(eigenvectors);

	std::cout << eigenvalues.transpose() << std::endl;

	return 0;
}