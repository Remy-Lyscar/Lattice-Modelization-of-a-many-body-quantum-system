#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>


#include "hamiltonian.h"
#include "operator.h"
#include "neighbours.h"



int main() {

	// PARAMETERS OF THE MODEL
	int m = 10; // number of sites
	int n = 10; // number of bosons 
	double J = 10; // hopping parameter
	double U = 0.1; // on-site interaction
	double mu = 0.1; // chemical potential
	Neighbours neighbours(m);
	neighbours.chain_neighbours(); // list of neighbours
	const std::vector<std::vector<int>>& nei = neighbours.getNeighbours();
	
	// INITIALIZATION OF THE HAMILTONIAN
	BH hmatrix (nei, m, n, J, U, mu);
	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();

	Operator H(std::move(smatrix));

	// DIAGONALIZATION OF THE HAMILTONIAN
	Eigen::VectorXd v_0 = Eigen::VectorXd::Random(H.size()); // Random initial vector
	Eigen::MatrixXd eigenvectors1;
	Eigen::MatrixXd eigenvectors3;
	Eigen::MatrixXd V;
	int k = 10; // Number of eigenvalues to calculate

	std::cout << "Dimension of the Hilbert space : " <<  H.size() << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::VectorXd eigenvalues1 = H.FOLM_eigen(k, v_0,eigenvectors1); // FOLM
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "FOLM execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "smallest eigenvalue : " << eigenvalues1.transpose()[0] << std::endl;
    std::cout << "number of calculated eigenvalues : " << eigenvalues1.size() << std::endl << std::endl;
	
	start = std::chrono::high_resolution_clock::now();
	Eigen::VectorXcd eigenvalues2 = H.IRLM_eigen(1); // IRLM
	end = std::chrono::high_resolution_clock::now();
	duration = end - start;
	std::cout << "IRLM execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "smallest eigenvalue : " << eigenvalues2.transpose() << std::endl;
    std::cout << "number of calculated eigenvalues : " << eigenvalues2.size() << std::endl << std::endl;

	start = std::chrono::high_resolution_clock::now();
	Eigen::VectorXd eigenvalues3 = H.exact_eigen(eigenvectors3); // Exact diagonalization
	end = std::chrono::high_resolution_clock::now();
	duration = end - start;
	std::cout << "Exact diagonalization execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "smallest eigenvalue : " << eigenvalues3.transpose()[0] << std::endl;
    std::cout << "number of calculated eigenvalues : " << eigenvalues3.size() << std::endl << std::endl;

	

	// PLOTTING THE EIGENVALUES
	// std::ofstream file("eigenvalues.txt");
	// if (file.is_open()) {
	// 	for (int i = 0; i < eigenvalues1.size(); ++i) {
	// 		file << eigenvalues1[i] << std::endl;
	// 	}
	// 	file.close();
	// }
	// else {
	// 	std::cerr << "Error when loading the file." << std::endl;
	// 	return 1;
	// }
	// // int result = system("python plot.py");
	// // if (result != 0) {
	// // 	std::cerr << "Error when executing Python script." << std::endl;
	// // 	return 1;
	// // }

	// PHASE TRANSITION CALCULATIONS
	double boson_density = H.boson_density(0.1);
	std::cout << "boson density : " << boson_density << std::endl;

	double compressibility = H.compressibility(0.1);
	std::cout << "compressibility : " << compressibility << std::endl;

	return 0;
}