#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
// #include <chrono>
// #include <complex>

#include "hamiltonian.h"
#include "operator.h"
#include "neighbours.h"



int main() {

	/// PARAMETERS OF THE MODEL

	int m = 7; // number of sites
	int n = 7; // number of bosons 
	double J = 100; // hopping parameter
	// double U = 0.1; // on-site interaction
	// double mu = 1; // chemical potential

	/// GEOMETRY OF THE LATTICE
	
	Neighbours neighbours(m);
	neighbours.chain_neighbours(); // list of neighbours
	const std::vector<std::vector<int>>& nei = neighbours.getNeighbours();
	
	/// INITIALIZATION OF THE HAMILTONIAN
	std::ofstream file("phase.txt");
	for (double mu=1; mu<102; mu+=10){
		for (double U=1; U<102; U+=10){
			BH hmatrix (nei, m, n, J, U, mu);
			Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
			Operator H(std::move(smatrix));
			double gap_ratio = H.gap_ratio();
			double boson_density = H.boson_density(1, n);
			double compressibility = H.compressibility(1, n);
			file << mu << " " << U << " " << gap_ratio << " " << boson_density << " " << compressibility << std::endl;
		}
	}
	file.close();

	
	/// DIAGONALIZATION
	
	// std::cout << "Dimension of the Hilbert space : " <<  H.size() << std::endl << std::endl;

	/// DIAGONALIZATION OF THE HAMILTONIAN
	
	// USING THE FOLM 
	// int k = 10; // Number of eigenvalues to calculate
	// Eigen::VectorXd v_0 = Eigen::VectorXd::Random(H.size()); // Random initial vector
	// Eigen::MatrixXd eigenvectors1;
	// Eigen::MatrixXd V;
	// auto start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXd eigenvalues1 = H.FOLM_eigen(k, v_0,eigenvectors1); // FOLM
	// auto end = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> duration = end - start;
	// std::cout << "FOLM execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << eigenvalues1.transpose()[0] << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues1.size() << std::endl << std::endl;
	
	// USING THE IRLM 
	// start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXcd eigenvalues2 = H.IRLM_eigen(1); // IRLM
	// end = std::chrono::high_resolution_clock::now();
	// duration = end - start;
	// std::cout << "IRLM execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << std::real(eigenvalues2.transpose()[0]) << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues2.size() << std::endl << std::endl;

	// // USING EXACT DIAGONALIZATION
	// Eigen::MatrixXd eigenvectors3;
	// start = std::chrono::high_resolution_clock::now();
	// Eigen::VectorXd eigenvalues3 = H.exact_eigen(eigenvectors3); // Exact diagonalization
	// end = std::chrono::high_resolution_clock::now();
	// duration = end - start;
	// std::cout << "Exact diagonalization execution time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "smallest eigenvalue : " << eigenvalues3.transpose()[0] << std::endl;
    // std::cout << "number of calculated eigenvalues : " << eigenvalues3.size() << std::endl << std::endl;


	/// PLOTTING THE EIGENVALUES
	// std::ofstream file("ground.txt");
	// if (file.is_open()) {
	// 	for (int i = 0; i < eigenvalues2.size(); ++i) {
	// 		file << eigenvalues2[i] << std::endl;
	// 	}
	// 	file.close();
	// }
	// else {
	// 	std::cerr << "Error when loading the file." << std::endl;
	// 	return 1;
	// }

	// int result = system("python3 plot.py");
	// if (result != 0) {
	// 	std::cerr << "Error when executing Python script." << std::endl;
	// 	return 1;
	// }


	/// PHASE TRANSITION CALCULATIONS
	// double boson_density = H.boson_density(0.1, n);
	// std::cout << "boson density : " << boson_density << std::endl;

	// double compressibility = H.compressibility(0.1, n);
	// std::cout << "isothermal compressibility : " << compressibility << std::endl << std::endl;

	return 0;
}