#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>


#include "hamiltonian.h"
#include "operator.h"
#include "neighbours.h"



int main() {

	// PARAMETERS OF THE MODEL
	int m = 3; // number of sites
	int n = 3; // number of bosons 
	double J = 1; // hopping parameter
	double U = 0; // on-site interaction
	double mu = 1; // chemical potential
	Neighbours neighbours(m);
	neighbours.chain_neighbours(m); // list of neighbours
	const std::vector<std::vector<int>>& nei = neighbours.getNeighbours();

	// INITIALIZATION OF THE HAMILTONIAN
	BH hmatrix (nei, m, n, J, U, mu);
	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
	Operator H(std::move(smatrix));


	// DIAGONALIZATION OF THE HAMILTONIAN
	Eigen::VectorXd v_0 = Eigen::VectorXd::Random(H.size()); // Random initial vector
	Eigen::MatrixXd eigenvectors1;
	Eigen::MatrixXd eigenvectors2;
	Eigen::MatrixXd V;
	int k = 10; // Number of eigenvalues to calculate

	Eigen::VectorXd eigenvalues1 = H.FOLM_eigen(k, v_0,eigenvectors1); // FOLM
	std::cout << eigenvalues1.transpose() << std::endl << eigenvalues1.size() << std::endl;
	 
	Eigen::VectorXd eigenvalues2 = H.exact_eigen(eigenvectors2); // Exact diagonalization
	std::cout << eigenvalues2.transpose() << std::endl << eigenvalues2.size() << std::endl;

	Eigen::VectorXcd eigenvalues3 = H.IRLM_eigen(k); // IRLM
	std::cout << eigenvalues3.transpose() << std::endl << eigenvalues3.size() << std::endl;


	// PLOTTING THE EIGENVALUES
	std::ofstream file("eigenvalues.txt");
	if (file.is_open()) {
		for (int i = 0; i < eigenvalues1.size(); ++i) {
			file << eigenvalues1[i] << std::endl;
		}
		file.close();
	}
	else {
		std::cerr << "Error when loading the file." << std::endl;
		return 1;
	}
	// int result = system("python plot.py");
	// if (result != 0) {
	// 	std::cerr << "Error when executing Python script." << std::endl;
	// 	return 1;
	// }

	return 0;
}