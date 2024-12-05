#include <Eigen/Dense>
#include <cmath>
#include <vector>

int sum(Eigen::VectorXd state, int index1, int index2) { // calculate the sum of the elements of a vector between 2 index
	int s = 0;
	for (int i = index1; i < index2 + 1; i++) {
		s += state[i];
	}
	return s;
}

bool next_lexicographic(Eigen::VectorXd state, int m, int n) { // calcule le prochain vecteur de la base
	if (state[m-1] == n) {
		return false;
	}
	for (int k = m-2; k > -1 ; k--) {
		if (state[k] != 0) {
			state[k] -= 1;
			state[k + 1] = n - sum(state, 0, k);
		}
	}
	return true;
}

Eigen::MatrixXd init_lexicographic(int m, int n) { // prend une matrice de taille adaptée pour y mettre les vecteurs de la base en colonne
	Eigen::MatrixXd basis; // initialize the matrix that stores the vector
	basis.resize(m, 0); 

	std::vector<int> s(m, 0);
	s[0] = n;

	Eigen::VectorXd<int> state = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
	s.clear();
	
	do {
		basis.conservativeResize(basis.rows(), basis.cols()+1); // agrandit la matrice d'une colonne
		basis.col(basis.cols()-1) = state;  // ajoute le vecteur à la dernière colonne de la matrice basis
	} while (next_lexicographic(state, m, n)) // calcule le prochain vecteur de la base
}

std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

double calculate_tag(const Eigen::MatrixXd basis, const std::vector primes, int k){ // calculate the tag of the kth column of the matrix
	double tag = 0;
	for (int i = 0; i < basis.rows(), ; i++) {
		tag += basis.coeff(i, k) * log(primes[i]); //see the formula page 7
	}
	return tag;
}

Eigen::VectorXd calculate_tags(const Eigen::MatrixXd basis, const std::vector primes) { // calculate and store the tags of each state of the basis
	std::vector<int> t(basis.cols(), 0);
	for (int i = 0; i < basis.cols(); i++) {
		t[i] = calculate_tag(basis, primes, i);
	}
	Eigen::VectorXd<int> tags = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(t.data(), t.size());
	t.clear();
	return tags;
}

voidsort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) { // sort by ascending order the tags and the vectors of the basis
	for (int i = 0; i < tags.size(); i++) {
		for (int j = i + 1; j < tags.size(); j++) {
			if (tags[j] < tags[i]) {
				std::swap(tags[i], tags[j]);
				basis.col(i).swap(basis.col(j));
			}
		}
	}
}

int search_tag(const Eigen::VectorXd tags, double x) { // gives the index of the wanted tag x by the Newton method
	int a = 0;
	int b = tags.size()-1;
	int m = (a + b) / 2;
	while (fabs(tags[m]-x) > 1e-3 and a <= b) {
		if (tags[m] < x) {
			a = m+1;
		}
		else {
			b = m-1;
		}
		mid = (a + b) / 2;
	}
	return m;
}

void fill_hopping(const Eigen::MatrixXd basis, const Eigen::VectorXd tags, const std::vector<std::vector<int>> neighbours, Eigen::SparseMatrixXd& hamiltonian, double J) { // CETTE METHODE N'EST PAS OPTIMISE
	for (int k = 0; k < basis.cols(); k++) {
		for (auto i : neighbours) {
			for (auto j : neighbours[i]) { // we want to calculate the term a_i^+ a_j and its complex conjugate
				Eigen::VectorXd state = basis.col(k).eval();
				if (basis.coeff(i, k) >= 1 and basis.coeff(j,k) >= 1) {
					state[i] += 1; 
					state[j] -= 1;
					double x = tag(state);
					int index = search_tag(tags, x);
					double value = sqrt( (basis.coeff(i, k)+1) * basis.coeff(j, k) );
					hamiltonian.insert(index, k) = -J * value;
					hamiltonian.insert(k, index) = -J * value;
				}
			}
		}
	}
}

void fill_interaction(const Eigen::MatrixXd basis, Eigen::SparseMatrixXd& hamiltonian, double U) { // CETTE METHODE N'EST PAS OPTIMISE
	for (int k = 0; k < basis.cols(); k++) {
		double value = 0;
		for (int i = 0; i < basis.rows(); i++) {
			double ni = basis.coeff(i, k);
			value += (ni + 1) * ni;
		}
		hamiltonian.insert(k, k) = U * value;
	}
}

void fill_chemical(const Eigen::MatrixXd basis, Eigen::SparseMatrixXd& hamiltonian, double mu) { // CETTE METHODE N'EST PAS OPTIMISE
	for (int k = 0; k < basis.cols(); k++) {
		double value = 0;
		for (int i = 0; i < basis.rows(); i++) {
			double ni = basis.coeff(i, k);
			value += ni;
		}
		hamiltonian.insert(k, k) = - mu * value;
	}
}

			