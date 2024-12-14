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
	int m = 4;
	int n = 4;

	double J = 1;
	double U = 0;
	double mu = 1;

	std::vector < std::vector<int>> neighbours = {{1,3}, {0,2}, {1,3},{2,0}};

	Hamiltonian hmatrix (neighbours, m, n, J, U, mu);

	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();

	Operator H(std::move(smatrix));
	Eigen::VectorXd v_0 = Eigen::VectorXd::Random(H.size());
	Eigen::MatrixXd eigenvectors;
	Eigen::MatrixXd V;
	Eigen::VectorXd eigenvalues = H.lanczos_eigen(10, v_0, V, eigenvectors);

	std::cout << eigenvalues.transpose() << std::endl;

	return 0;
}