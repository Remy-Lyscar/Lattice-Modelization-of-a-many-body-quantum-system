#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <vector>

#include "hamiltonian.h"


// ELEMENTARY FUNCTIONS

/* calculate the sum of the elements of a vector between 2 index */
int Hamiltonian::sum(const Eigen::VectorXd& state, int index1, int index2) const { 
	int s = 0;
	for (int i = index1; i <= index2; i++) {
		s += state[i];
	}
	return s;
}


// INITIALIZE THE HILBERT SPACE BASIS

/* calculate the next state of the Hilbert space in lexicographic order */
bool Hamiltonian::next_lexicographic(Eigen::VectorXd& state, int m, int n) const {
	for (int k = m - 2; k > -1; k--) {
		if (state[k] != 0) {
			state[k] -= 1;
			state[k + 1] = n - sum(state, 0, k);
			for (int i = k + 2; i < m; i++) {
				state[i] = 0;
			}
			return true;
		}
	}
	return false;
}

/* creates a matrix that has the vectors of the Hilbert space basis in columns */
Eigen::MatrixXd Hamiltonian::init_lexicographic(int m, int n) const {
	Eigen::MatrixXd basis(m, 1);
	basis.col(0).setZero();
	basis(0, 0) = n;

	Eigen::VectorXd state = basis.col(0);

	while (next_lexicographic(state, m, n)) {
		basis.conservativeResize(Eigen::NoChange, basis.cols() + 1);
		basis.col(basis.cols() - 1) = state;
	}
	return basis;
}


// SORT THE HILBERT SPACE BASIS TO FACILITATE CALCULUS

/* calculate the unique tag of the kth column of the matrix */
double Hamiltonian::calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k) const {
	double tag = 0;
	for (int i = 0; i < basis.rows(); i++) {
		tag += basis.coeff(i, k) * log(primes[i]);
	}
	return tag;
}

/* calculate and store the tags of each state of the Hilbert space basis */
Eigen::VectorXd Hamiltonian::calculate_tags(const Eigen::MatrixXd& basis, const std::vector<int>& primes) const {
	Eigen::VectorXd tags(basis.cols());
	for (int i = 0; i < basis.cols(); i++) {
		tags[i] = calculate_tag(basis, primes, i);
	}
	return tags;
}

/* sort the states of the Hilbert space by ascending order compared by their tags*/
void Hamiltonian::sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) const {
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
int Hamiltonian::search_tag(const Eigen::VectorXd& tags, double x) const {
	int a = 0;
	int b = tags.size() - 1;
	int m = (a + b) / 2;
	while (fabs(tags[m] - x) > 1e-3 && a <= b) {
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

/* fill the hopping term of the Hamiltonian */
void Hamiltonian::fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, const std::vector<int>& primes, Eigen::SparseMatrix<double>& hmatrix, double J) const {
	std::vector<Eigen::Triplet<double>> tripletList;
	for (int k = 0; k < basis.cols(); k++) {
		for (int i = 0; i < neighbours.size(); i++) {
			for (int j = 0; j < neighbours[i].size(); j++) {
				Eigen::VectorXd state = basis.col(k).eval();
				if (basis.coeff(i, k) >= 1 && basis.coeff(j, k) >= 1) {
					state[i] += 1;
					state[j] -= 1;
					double x = calculate_tag(state, primes, i);
					int index = search_tag(tags, x);
					double value = sqrt((basis.coeff(i, k) + 1) * basis.coeff(j, k));
					tripletList.push_back(Eigen::Triplet<double>(index, k, -J * value));
					tripletList.push_back(Eigen::Triplet<double>(k, index, -J * value));
				}
			}
		}
	}
	hmatrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

/* fill the interaction term of the Hamiltonian */
void Hamiltonian::fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double U) const {
	std::vector<Eigen::Triplet<double>> tripletList;

	for (int k = 0; k < hmatrix.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(hmatrix, k); it; ++it) {
			tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < basis.cols(); k++) {
		double value = 0;
		for (int i = 0; i < basis.rows(); i++) {
			double ni = basis.coeff(i, k);
			value += (ni + 1) * ni;
		}
		tripletList.push_back(Eigen::Triplet<double>(k, k, U * value));
	}
	hmatrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

/* fill the chemical potential term of the Hamiltonian */
void Hamiltonian::fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double mu) const {
	std::vector<Eigen::Triplet<double>> tripletList;
	
	for (int k = 0; k < hmatrix.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(hmatrix, k); it; ++it) {
			tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < basis.cols(); k++) {
		double value = 0;
		for (int i = 0; i < basis.rows(); i++) {
			double ni = basis.coeff(i, k);
			value += ni;
		}
		tripletList.push_back(Eigen::Triplet<double>(k, k, -mu * value));
	}
	hmatrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

// CONSTRUCTOR

/* constructor by default for the Hamiltonian class */
Hamiltonian::Hamiltonian(const std::vector<std::vector<int>>& neighbours, int m, int n, double J, double U, double mu) : neighbours(neighbours), m(m), n(n), J(J), U(U), mu(mu) {
    Eigen::MatrixXd basis = init_lexicographic(m, n);
    mat.resize(basis.cols(), basis.cols());
    mat.setZero();
    if (J != 0) {
        std::vector<int> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
        Eigen::VectorXd tags = calculate_tags(basis, primes);
        sort_basis(tags, basis);
        fill_hopping(basis, tags, neighbours, primes, mat, J);
    }
    if (U != 0) {
        fill_interaction(basis, mat, U);
    }
    if (mu != 0) {
        fill_chemical(basis, mat, mu);
    }
}

//UTILITY FUNCTIONS

/* return the Hamiltonian sparse matrix */
Eigen::SparseMatrix<double> Hamiltonian::getHamiltonian() const {
	return mat;
}