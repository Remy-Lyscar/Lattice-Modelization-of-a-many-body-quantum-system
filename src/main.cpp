#include <pybind11/pybind11.h>  // includes also Python.h -> it has to be the first include!
                                // Indeed, Python.h defines some preprocessor variables
                                // that may affect the behavior of the standard headers.
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>  //python interpreter
#include <pybind11/stl.h>  // type conversion

/* Pybind11 is a library that exposes C++ type in Python and vice versa in order to make bindings.
It is lighter than the Boost.Python project (that was designed for the same purpose),
because it only matches Cpp11 compilers and above.*/

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



namespace py = pybind11;


int main(){

	std::cout << "/*----- Starting pybind -----*/" << std::endl;
	py::scoped_interpreter guard{}; // start python interpreter, dies when out of scope

	auto math = py::module::import("math");
    double root_two = math.attr("sqrt")(2.0).cast<double>();

    std::cout << "The square root of 2 is: " << root_two << "\n";


	py::object eigsh = py::module::import("scipy.sparse.linalg").attr("eigsh");


	int N = 2;  

	XY XY(3,3);


	// auto result = eigsh(M, py::arg("k") = 1, py::arg("which") = "SA", py::arg("return_eigenvectors") = false);
	// auto vaps = result.cast<std::array<std::complex<double>, 1>>();

	return 0; 
}



// int main() {

// 	// PARAMETERS OF THE MODEL
// 	int m = 4; // number of sites
// 	int n = 4; // number of bosons 
// 	double J = 1; // hopping parameter
// 	double U = 0; // on-site interaction
// 	double mu = 1; // chemical potential
// 	std::vector < std::vector<int>> neighbours = Hamiltonian::chain_neighbours(m); // list of neighbours


// 	// INITIALIZATION OF THE HAMILTONIAN
// 	Hamiltonian hmatrix (neighbours, m, n, J, U, mu);
// 	Eigen::SparseMatrix<double> smatrix = hmatrix.getHamiltonian();
// 	Operator H(std::move(smatrix));
// 	std::cout << H.size() << std::endl;

// 	// INITIALIZATION OF THE RANDOM VECTOR FOR LANCZOS ALGORITHM
// 	RandomVector randomVector(H.size(), n, 42);
// 	Eigen::VectorXd v_0 = randomVector.uniform();


// 	// DIAGONALIZATION OF THE HAMILTONIAN
// 	Eigen::MatrixXd eigenvectors1;
// 	Eigen::MatrixXd eigenvectors2;
// 	Eigen::MatrixXd V;

// 	Eigen::VectorXd eigenvalues1 = H.lanczos_eigen(H.size(), v_0, V, eigenvectors1); // Lanczos diagonalization
// 	std::cout << eigenvalues1.transpose() << std::endl << eigenvalues1.size() << std::endl;
	 
// 	Eigen::VectorXd eigenvalues2 = H.exact_eigen(eigenvectors2); // Exact diagonalization
// 	std::cout << eigenvalues2.transpose() << std::endl << eigenvalues2.size() << std::endl;

// 	// PLOTTING THE EIGENVALUES
// 	std::ofstream file("eigenvalues.txt");
// 	if (file.is_open()) {
// 		for (int i = 0; i < eigenvalues1.size(); ++i) {
// 			file << eigenvalues1[i] << std::endl;
// 		}
// 		file.close();
// 	}
// 	else {
// 		std::cerr << "Error when loading the file." << std::endl;
// 		return 1;
// 	}
// 	int result = system("python plot.py");
// 	if (result != 0) {
// 		std::cerr << "Error when executing Python script." << std::endl;
// 		return 1;
// 	}

// 	return 0;
// }