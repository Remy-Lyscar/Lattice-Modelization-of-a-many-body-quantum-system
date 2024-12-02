#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <iostream>
#include <stdexcept>

class Hamiltonian {
private:

	Eigen::SparseMatrix<double> smat;
	
	double ref; // threshold under which a value is considered null

	/* calculate the dimension of the Hilbert space */
	long long m,n,D; // number of sites, bosons and Hilbert space dimension;

	long long factorial(long long n) { // calculate factorial n with a recursive function
		if (n == 0) {
			return 1;
		}
		else {
			return n * factorial(n - 1);
		}
	}
	long long dimension(long long m, long long n) {// calculate the dimension of the Hilbert space for n bosons on m sites
		return factorial(n + m + 1) / (factorial(n) * factorial(m - 1));
	}
	

	/* implement the Lanczos algorithm for a sparse matrix for nb_iter iterations starting with vector v_0 */
	void lanczos_diag(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& T, Eigen::MatrixXd& V) {

		T.resize(nb_iter, nb_iter); //resize the matrix T to a matching size
		V.resize(D, nb_iter); //resize the matrix V to a matching size

		v_0.normalize(); // normalize the starting vector v_0
		V.col(0) = v_0; // put the first vector v_0 in the matrix V

		double alpha, beta;

		//Lanczos algorithm for nb_iter iterations
		for (int i = 0; i < nb_iter; i++) {
			Eigen::VectorXd w = this->smat * V.col(i); // calculate the next vector w 
			alpha = (V.col(i)).dot(w);
			for (int j = 0; j < i; j++) {
				w = w - (V.col(j)).dot(w) * V.col(j); // orthogonalize the vector w with respect to the previous vectors of the Krylov basis
			}
			beta = w.norm(); // calculate the norm of the vector w
			if (beta < ref) {
				break; // if beta is null or almost null the algorithm stops
			}
			else {
				w = w / beta; // normalize the vector w
				V.col(i + 1) = w; // add the vector w to the matrix V of vectors of the Krylov basis
				T(i, i) = alpha; // add the ith diagonal element of the tridiagonal matrix T
				if (i > 0) {
					T(i, i - 1) = beta; // add the ith non-diagonal element of the tridiagonal matrix T
					T(i - 1, i) = beta; // add the ith non-diagonal element of the tridiagonal matrix T
				}
			}
		}
	}
	
	void sort_eigen(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) { // A CHANGER
		std::vector<int> indices(eigenvalues.size());
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(), [&eigenvalues](int i, int j) {
			return eigenvalues[i] > eigenvalues[j];
			});

		Eigen::VectorXd sortedEigenvalues(eigenvalues.size());
		Eigen::MatrixXd sortedEigenvectors(eigenvectors.rows(), eigenvectors.cols());

		for (int i = 0; i < indices.size(); ++i) {
			sortedEigenvalues[i] = eigenvalues[indices[i]];
			sortedEigenvectors.col(i) = eigenvectors.col(indices[i]);
		}

		eigenvalues = sortedEigenvalues;
		eigenvectors = sortedEigenvectors;
	}

public:
	Hamiltonian(long long m, long long n, double threshold = 1e-6) : m(m), n(n), ref(threshold), D(dimension(m, n)), smat(D, D) {}; // A CHANGER

	/* add a matrix to an operand of type SparseMatrix with same size */
	Hamiltonian operator + (const Hamiltonian& operand) {
		if (this->smat.rows() != operand.smat.rows() or this->smat.cols() != operand.smat.cols()) { // verify that the operands have matching size
			throw std::invalid_argument("Matrix should have matching size.");
		}
		Hamiltonian result(operand.smat.rows(), operand.smat.cols, ref); //A CHANGER
		result.smat = (this->smat + operand.smat).pruned(ref);
		return result; // removes element smaller than ref
	}

	/* multiply a sparse matrix by a multiplicand of type SparseMatrix with same size */
	Hamiltonian operator * (const Hamiltonian& multiplicand) {
		if (this->smat.cols() != multiplicand.smat.rows()) {
			throw std::invalid_argument("Number of columns of multiplier must equal number of rows of multiplicand."); // verify that the operands have matching size
		}
		Hamiltonian result(this->smat.rows(), multiplicand.smat.cols(), ref); //A CHANGER
		result.smat = (this->smat * multiplicand.smat).pruned(ref);
		return result; // removes element smaller than ref
	}

	/* multiply a sparse matrix by a vector with concomitant size */
	Eigen::VectorXd operator * (const Eigen::VectorXd& vector) {
		if (this->smat.rows() != vector.size()) { // verify that the operands have matching size
			throw std::invalid_argument("Number of matrix rows must equal vector size.");
		}
		return this->smat * vector;
	}

	/* calculate the approximate eigenvalues and eigenvectors of the hamiltonian in Krylov space using the Lanczos algorithm */
	Eigen::VectorXd lanczos_eigen(int nb_iter, Eigen::VectorXd& v_0, Eigen::MatrixXd& V) {
		Eigen::MatrixXd T; // initialize the tridiagonal matrix T 
		lanczos_diag(nb_iter, v_0, T, V); // tridiagonalize the hamiltonian using Lanczos algorithm

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(T); // solve the eigen problem for T
		if (eigensolver.info() != Eigen::Success) { // verify if the eigen search is a success
			throw std::runtime_error("Eigenvalue computation failed.");
		}
		Eigen::MatrixXd eigenvectors = V * eigensolver.eigenvectors(); // eigenvectors of the hamiltonian
		Eigen::VectorXd eigenvalues = eigensolver.eigenvalues(); // eigenvalues of the hamiltonian
		sort_eigen(eigenvalues, eigenvectors);
		return eigenvalues;
	}

	/* calculate the exact eigenvalues and eigenvectors of the hamiltonian by an exact diagonalization */
	Eigen::VectorXd exact_eigen(Eigen::MatrixXd& eigenvectors) {
		Eigen::SparseSelfAdjointEigenSolver<Eigen::SparseMatrix> eigensolver(this->smat); // solve the eigen problem for the hamiltonian
		if (eigensolver.info() != Eigen::Success) { // verify if the eigen search is a success
			throw std::runtime_error("Eigenvalue computation failed.");
		}
		eigenvectors = eigensolver.eigenvectors(); // eigenvectors of the hamiltonian
		Eigen::VectorXd eigenvalues = eigensolver.eigenvalues(); // eigenvalues of the hamiltonian
		sort_eigen(eigenvalues, eigenvectors);
		return eigenvalues;
	}
}






	
	




