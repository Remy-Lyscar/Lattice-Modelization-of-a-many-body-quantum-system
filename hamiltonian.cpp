#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <vector>

//#include "hamiltonian.h"

class Hamiltonian {
private:

	Eigen::SparseMatrix<double> mat;
	std::vector<std::vector<int>> neighbours;
	int m, n;
	double J, U, mu;


	//  DIMENSION OF THE HILBERT SPACE

	int factorial(int n) { // calculate factorial n with a recursive function
		if (n == 0) {
			return 1;
		}
		else {
			return n * factorial(n - 1);
		}
	}

	int dimension(int m, int n) {// calculate the dimension of the Hilbert space for n bosons on m sites
		return factorial(n + m + 1) / (factorial(n) * factorial(m - 1));
	}


	// ELEMENTARY FUNCTIONS

	int sum(const Eigen::VectorXd& state, int index1, int index2) { // calculate the sum of the elements of a vector between 2 index
		int s = 0;
		for (int i = index1; i < index2 + 1; i++) {
			s += state[i];
		}
		return s;
	}


	// INITIALIZE THE HILBERT SPACE BASIS

	/* calculate the next state of the Hilbert space in lexicographic order */
	bool next_lexicographic(Eigen::VectorXd& state, int m, int n) {
		if (state[m - 1] == n) {
			return false;
		}
		for (int k = m - 2; k > -1; k--) {
			if (state[k] != 0) {
				state[k] -= 1;
				state[k + 1] = n - sum(state, 0, k);
			}
		}
		return true;
	}

	/* creates a matrix that has the vectors of the Hilbert space basis in columns */
	Eigen::MatrixXd init_lexicographic(int m, int n) {
		Eigen::MatrixXd<double> basis; // initialize the matrix that will store the states
		basis.resize(m, 0);

		std::vector<int> s(m, 0);
		s[0] = n;

		Eigen::VectorXd state = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
		s.clear();

		do {
			basis.conservativeResize(basis.rows(), basis.cols() + 1); // add one column to the matrix
			basis.col(basis.cols() - 1) = state;  // add the last calculated state to the matrix
		} while (next_lexicographic(state, m, n)); // calculate the next state of the basis
		return basis;
	}


	// SORT THE HILBERT SPACE BASIS TO FACILITATE CALCULUS

	/* calculate the unique tag of the kth column of the matrix */
	double calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k) {
		double tag = 0;
		for (int i = 0; i < basis.rows(); i++) {
			tag += basis.coeff(i, k) * log(primes[i]); //see the formula page 7
		}
		return tag;
	}

	/* calculate and store the tags of each state of the Hilbert space basis */
	Eigen::VectorXd calculate_tags(const Eigen::MatrixXd& basis, const std::vector& primes) {
		std::vector<int> t(basis.cols(), 0);
		for (int i = 0; i < basis.cols(); i++) {
			t[i] = calculate_tag(basis, primes, i);
		}
		Eigen::VectorXd tags = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(t.data(), t.size());
		t.clear();
		return tags;
	}

	/* sort the states of the Hilbert space by ascending order compared by their tags*/
	void sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) {
		for (int i = 0; i < tags.size(); i++) {
			for (int j = i + 1; j < tags.size(); j++) {
				if (tags[j] < tags[i]) {
					std::swap(tags[i], tags[j]);
					basis.col(i).swap(basis.col(j));
				}
			}
		}
	}


	// FILL THE HAMILTONIAN OF THE SYSTEM

	/* gives the index of the wanted tag x by the Newton method */
	int search_tag(const Eigen::VectorXd& tags, double x) {
		int a = 0;
		int b = tags.size() - 1;
		int m = (a + b) / 2;
		while (fabs(tags[m] - x) > 1e-3 and a <= b) {
			if (tags[m] < x) {
				a = m + 1;
			}
			else {
				b = m - 1;
			}
			m = (a + b) / 2;
		}
		return m;
	}

	void fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, Eigen::SparseMatrixXd<double>& hmatrix, double J) { // CETTE METHODE N'EST PAS OPTIMISE
		for (int k = 0; k < basis.cols(); k++) {
			for (int i = 0; i < neighbours.size(); i++) {
				for (int j = 0; j < neighbours[i].size(); j++) { // we want to calculate the term a_i^+ a_j and its complex conjugate
					Eigen::VectorXd state = basis.col(k).eval();
					if (basis.coeff(i, k) >= 1 and basis.coeff(j, k) >= 1) {
						state[i] += 1;
						state[j] -= 1;
						double x = tag(state);
						int index = search_tag(tags, x);
						double value = sqrt((basis.coeff(i, k) + 1) * basis.coeff(j, k));
						hmatrix.insert(index, k) = -J * value;
						hmatrix.insert(k, index) = -J * value;
					}
				}
			}
		}
	}

	void fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrixXd<double>& hmatrix, double U) { // CETTE METHODE N'EST PAS OPTIMISE
		for (int k = 0; k < basis.cols(); k++) {
			double value = 0;
			for (int i = 0; i < basis.rows(); i++) {
				double ni = basis.coeff(i, k);
				value += (ni + 1) * ni;
			}
			hmatrix.insert(k, k) = U * value;
		}
	}

	void fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrixXd<double>& hmatrix, double mu) { // CETTE METHODE N'EST PAS OPTIMISE
		for (int k = 0; k < basis.cols(); k++) {
			double value = 0;
			for (int i = 0; i < basis.rows(); i++) {
				double ni = basis.coeff(i, k);
				value += ni;
			}
			hmatrix.insert(k, k) = -mu * value;
		}
	}

public:
	Hamiltonian(const std::vector<std::vector<int>>& neighbours, int m, int n, double J, double U, double mu) : neighbours(neighbours), m(m), n(n), J(J), U(U), mu(mu) {
		Eigen::MatrixXd basis = init_lexicographic(m, n);
		mat.resize(basis.cols(), basis.cols());
		mat.setZero();
		if (J != 0) {
			std::vector<int> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
			Eigen::VectorXd tags = calculate_tags(basis, primes);
			sort_basis(tags, basis);
			fill_hopping(basis, tags, neighbours, mat, J);

		}
		if (U != 0) {
			fill_interaction(basis, mat, U);
		}
		if (mu != 0) {
			fill_chemical(basis, mat, mu);
		}
	}

	Eigen::SparseMatrix<double> getHamiltonian() {
		return mat;
	}
};