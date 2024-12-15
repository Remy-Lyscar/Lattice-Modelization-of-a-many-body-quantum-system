#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "hamiltonian.h"
#include "operator.h"
#include "randomvector.h"

/* generate the list of neighbours for a 1D chain */
std::vector<std::vector<int>> chain_neighbours(int m, int closed = 1) { // closed = 0 for open boundary conditions, closed = 1 for periodic boundary conditions
	std::vector<std::vector<int>> neighbours(m);
	for (int i = 0; i < m; ++i) {
		if (i > 0) {
			neighbours[i].push_back(i - 1); // Left neighbour
		}
		if (i < m - 1) {
			neighbours[i].push_back(i + 1); // Right neighbour
		}
	}
	if (closed == 1) { // Periodic boundary conditions
		neighbours[0].push_back(m - 1); 
		neighbours[m - 1].push_back(0); 
	}
	return neighbours;
}

int main() {

	// PARAMETERS OF THE MODEL
	int m = 4; // number of sites
	int n = 4; // number of bosons 
	double J = 1; // hopping parameter
	double U = 0; // on-site interaction
	double mu = 1; // chemical potential
	std::vector < std::vector<int>> neighbours = chain_neighbours(m); // list of neighbours


	// INITIALIZATION OF THE HAMILTONIAN
	Hamiltonian hmatrix (neighbours, m, n, J, U, mu);
	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
	Operator H(std::move(smatrix));
	std::cout << H.size() << std::endl;

	// INITIALIZATION OF THE RANDOM VECTOR FOR LANCZOS ALGORITHM
	RandomVector randomVector(H.size(), n, 42);
	Eigen::VectorXd v_0 = randomVector.uniform();


	// DIAGONALIZATION OF THE HAMILTONIAN
	Eigen::MatrixXd eigenvectors1;
	Eigen::MatrixXd eigenvectors2;
	Eigen::MatrixXd V;

	Eigen::VectorXd eigenvalues1 = H.lanczos_eigen(H.size(), v_0, V, eigenvectors1); // Lanczos diagonalization
	std::cout << eigenvalues1.transpose() << std::endl << eigenvalues1.size() << std::endl;
	 
	Eigen::VectorXd eigenvalues2 = H.exact_eigen(eigenvectors2); // Exact diagonalization
	std::cout << eigenvalues2.transpose() << std::endl << eigenvalues2.size() << std::endl;

	// PLOTTING THE EIGENVALUES
	std::ofstream file("eigenvalues.txt");
	if (file.is_open()) {
		for (int i = 0; i < eigenvalues.size(); ++i) {
			file << eigenvalues[i] << std::endl;
		}
		file.close();
	}
	else {
		std::cerr << "Error when loading the file." << std::endl;
		return 1;
	}
	int result = system("python plot.py");
	if (result != 0) {
		std::cerr << "Error when executing Python script." << std::endl;
		return 1;
	}

	return 0;
}