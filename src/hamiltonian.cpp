#include<iostream>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/SparseCore>
#include<complex>
#include<vector>

#include "hamiltonian.h"




///// IMPLEMENTATION OF HAMILTONIAN CLASS METHODS /////


    /// NEIGHBOURS ///

/* generate the list of neighbours for a 1D chain */
std::vector<std::vector<int>> Hamiltonian::chain_neighbours(int m, bool closed) { // by default closed = true for periodic boundary conditions, closed = false for open boundary conditions
	std::vector<std::vector<int>> neighbours(m);
	for (int i = 0; i < m; ++i) {
		if (i > 0) {
			neighbours[i].push_back(i - 1); // Left neighbour
		}
		if (i < m - 1) {
			neighbours[i].push_back(i + 1); // Right neighbour
		}
	}
	if (closed) { // Periodic boundary conditions
		neighbours[0].push_back(m - 1); 
		neighbours[m - 1].push_back(0); 
	}
	return neighbours;
}

/* generate the list of neighbours for a 2D square lattice */
std::vector<std::vector<int>> Hamiltonian::square_neighbours(int m, bool closed) { // by default closed = true for periodic boundary conditions, closed = false for open boundary conditions
    int side = static_cast<int>(std::sqrt(m));
    if (side * side != m) {
        throw std::invalid_argument("The number of sites (m) must be a perfect square.");
    }
    std::vector<std::vector<int>> neighbours(m);
    for (int i = 0; i < side; ++i) {
        for (int j = 0; j < side; ++j) {
            int index = i * side + j;
            if (j > 0) {
                neighbours[index].push_back(index - 1); // Left neighbour
            } else if (closed) {
                neighbours[index].push_back(index + side - 1);
            }
            if (j < side - 1) {
                neighbours[index].push_back(index + 1); // Right neighbour
            } else if (closed) {
                neighbours[index].push_back(index - side + 1);
            }
            if (i > 0) {
                neighbours[index].push_back(index - side); // Top neighbour
            } else if (closed) {
                neighbours[index].push_back(index + (side - 1) * side);
            }
            if (i < side - 1) {
                neighbours[index].push_back(index + side); // Bottom neighbour
            } else if (closed) {
                neighbours[index].push_back(index - (side - 1) * side);
            }
        }
    }
    return neighbours;
}


    /// DISPLAY FUNCTIONS ///

/* Display a sparse matrix in a standard form in the terminal, as a full matrix */
template <typename T> void Hamiltonian::displaySparseMatrix(const Eigen::SparseMatrix<T>& M) const{
    std::cout << "Displaying the sparse matrix at memory adress: " << &M << "\n \n" << std::endl; 
    for(int i = 0; i < M.rows(); i++){
        for(int j = 0; j < M.cols(); j++){
            T c = M.coeff(i, j);
            if (std::norm(c) < epsilon){
                std::cout << "(0,0)" << " "; 
            }
            else{
                std::cout  << "\033[34m" << c << " "; // Here we print the non zero elements in blue 
                                                      // so they are easy to see in the terminal
                std::cout << "\033[0m"; // Reset the color to the default one
            } 
        }
        std::cout << std::endl; 
    }
    std::cout << "\n \n" << std::endl;
}

/* Explicit template instantiation for the types you need */
template void Hamiltonian::displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const;
template void Hamiltonian::displaySparseMatrix(const Eigen::SparseMatrix<double>& M) const;


/*Display a dense matrix in the terminal */
void Hamiltonian::displayMatrix(const Eigen::MatrixXd& M) const{
    std::cout << "Displaying the matrix at memory adress: " << &M << "\n \n" << std::endl; 
    for(int i = 0; i < M.rows(); i++){
        for(int j = 0; j < M.cols(); j++){
            std::cout << M(i, j) << " "; 
        }
        std::cout << std::endl; 
    }
    std::cout << "\n \n" << std::endl;
}





///// IMPLEMENTATION OF THE XY CLASS METHODS /////


    /// CONSTRUCTOR ///

/* Constructor for the XY model with N sites and coupling parameter J = 1.0 */
XY::XY(int N_): N(N_), D(1<<N_), J(1.0), mu(0.0),H(D,D), S_x(2,2), S_y(2,2), S_z(2,2), I(2,2){
    
    using namespace std::complex_literals; 

    // Identity matrix
    I.setIdentity(); 

    // Initialization of the Pauli Matrices
    S_x.insert(0, 0) = 0.; 
    S_x.insert(1,0) = 1.; 
    S_x.insert(0,1) = 1.; 
    S_x.insert(1,1) = 0.;
    S_x = 0.5*S_x; // Spin operator is 1/2 times the Pauli matrix

    S_y.insert(0, 0) = 0.; 
    S_y.insert(1,0) = 0. + 1i; 
    S_y.insert(0,1) = 0. - 1i; 
    S_y.insert(1,1) = 0.; 
    S_y = 0.5*S_y; // Spin operator is 1/2 times the Pauli matrix

    S_z.insert(0, 0) = 1.; 
    S_z.insert(1,0) = 0.; 
    S_z.insert(0,1) = 0.; 
    S_z.insert(1,1) = -1.;
    S_z = 0.5*S_z; // Spin operator is 1/2 times the Pauli matrix

    computeHamiltonianXY(); // Initialization of the Hamiltonian of the XY model
    displaySparseMatrix(H); // Display the Hamiltonian of the XY model
}


    /// Private methods ///

/* Computes the Kronecker Product of two Eigen::SparseMatrix objects and stores the result in a Eigen::SparseMAtrix object*/
Eigen::SparseMatrix<std::complex<double>> XY::kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::SparseMatrix<std::complex<double>>& B) const{
    typedef Eigen::SparseMatrix<std::complex<double>> SparseMatrix;
    typedef Eigen::Triplet<std::complex<double>> Triplet; // The Triplet class is a utility class used to build SparseMatrix easily
    std::vector<Triplet> tripletList;
    tripletList.reserve(A.nonZeros() * B.nonZeros());

    // Can be simplified later if we only deal with square matrices
    int rowsA = A.rows();
    int colsA = A.cols();
    int rowsB = B.rows();
    int colsB = B.cols();

    SparseMatrix C(rowsA * rowsB, colsA * colsB);

    /* The algorithm is based on the idea: 
    for each pair of non zero elements of A and B, we compute the element that results from those in the C matrix
    Then the SparseMatrix C is created with the method setFromTriplets */

    for (int k = 0; k < A.outerSize(); ++k) {
        for (typename SparseMatrix::InnerIterator itA(A, k); itA; ++itA) {
            for (int l = 0; l < B.outerSize(); ++l) {
                for (typename SparseMatrix::InnerIterator itB(B, l); itB; ++itB) {
                    tripletList.push_back(Triplet(
                        itA.row() * rowsB + itB.row(),
                        itA.col() * colsB + itB.col(),
                        itA.value() * itB.value()
                    ));
                }
            }
        }
    }
    C.setFromTriplets(tripletList.begin(), tripletList.end());
    return C;
}

/* Method that computes the Kronecker Product of a Pauli matrix S at site i with the identity matrix at the other sites */
Eigen::SparseMatrix<std::complex<double>> XY::kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, int i) const{
    Eigen::SparseMatrix<std::complex<double>> result = (i==0) ? S: I;  
    for(int j = 1; j < N; j++){
        if(j == i){
            result = kroneckerProductSparse(result, S); 
        }
        else{
            result = kroneckerProductSparse(result, I); 
        }
    }
    return result; 
}

/* This method computes the full hamiltonian of the XY quantum system */
void XY::computeHamiltonianXY() {
    // Initialize H with the first term. Operator * is overloaded for SparseMatrix multiplication
    H = kroneckerPauli(S_x,0)*kroneckerPauli(S_x, 1) + kroneckerPauli(S_y,0)*kroneckerPauli(S_y, 1); 
    for(int i = 1; i < N-1; i++){
        H += kroneckerPauli(S_x,i)*kroneckerPauli(S_x, i+1) + kroneckerPauli(S_y,i)*kroneckerPauli(S_y, i+1);
    }
}


    /// Public methods ///

/* Return the Hamiltonian sparse matrix */
Eigen::SparseMatrix<std::complex<double>> XY::getHamiltonian() const {
	return H;
}

/* Display all the attributes of the XY model */
void XY::display_all() const {
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 
    displaySparseMatrix(H);
}





/////  IMPLEMENTATION OF THE BH CLASS METHODS  /////

    
    /// ELEMENTARY FUNCTIONS ///

/* Calculate the sum of the elements of a vector between 2 index */
int BH::sum(const Eigen::VectorXd& state, int index1, int index2) const { 
	int s = 0;
	for (int i = index1; i <= index2; i++) {
		s += state[i];
	}
	return s;
}


    /// DIMENSION OF THE HILBERT SPACE ///

/* Calculate the binomial coefficient */
int BH::binomial(int n, int k) const{
	if (k==0 || k==n){
		return 1;
	}
	if (k > n/2) {
		return binomial(n,n-k);
	}
	else{
		return n*binomial(n-1,k-1)/k;
	}
}

/* Calculate the dimension of the Hilbert space for n bosons on m sites */
int BH::dimension(int m, int n) const{
	return binomial(m + n - 1, m);
}


    /// INITIALIZE THE HILBERT SPACE BASIS ///

/* Calculate the next Fock state of the Hilbert space in lexicographic order */
bool BH::next_lexicographic(Eigen::VectorXd& state, int m, int n) const {
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

/* Create the matrix that has the Fock states of the Hilbert space basis in columns */
Eigen::MatrixXd BH::init_lexicographic(int m, int n) const {
	Eigen::MatrixXd basis(m, D);
	Eigen::VectorXd state = Eigen::VectorXd::Zero(m);
    state(0) = n;
	int col = 0;
	do {
		basis.col(col++) = state;
	} while (next_lexicographic(state, m, n)) ;
	return basis;
}


    /// SORT THE HILBERT SPACE BASIS TO FACILITATE CALCULUS ///

/* Calculate the unique tag of the kth column of the matrix */
double BH::calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k) const {
	double tag = 0;
	for (int i = 0; i < basis.rows(); i++) {
		tag += basis.coeff(i, k) * log(primes[i]);
	}
	return tag;
}

/* Calculate and store the tags of each state of the Hilbert space basis */
Eigen::VectorXd BH::calculate_tags(const Eigen::MatrixXd& basis, const std::vector<int>& primes) const {
	Eigen::VectorXd tags(basis.cols());
	for (int i = 0; i < basis.cols(); i++) {
		tags[i] = calculate_tag(basis, primes, i);
	}
	return tags;
}

/* Sort the states of the Hilbert space by ascending order compared by their tags*/
void BH::sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) const {
    std::vector<int> indices(tags.size());
    for (int i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&tags](int a, int b) {return tags[a] < tags[b];});
    for (int i = 0; i < indices.size(); ++i) {
        while (indices[i] != i) {
            int j = indices[i];
            std::swap(tags[i], tags[j]);
            basis.col(i).swap(basis.col(j));
            std::swap(indices[i], indices[j]);
        }
    }
}

/* Gives the index of the wanted tag x by the Newton method */
int BH::search_tag(const Eigen::VectorXd& tags, double x) const {
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


    /// FILL THE HAMILTONIAN OF THE SYSTEM ///

/* Fill the hopping term of the Hamiltonian */
void BH::fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, const std::vector<int>& primes, Eigen::SparseMatrix<double>& hmatrix, double J) const {
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(basis.cols() * basis.rows() * neighbours.size());
	for (int k = 0; k < basis.cols(); k++) {
		for (int i = 0; i < neighbours.size(); i++) {
			for (int j = 0; j < neighbours[i].size(); j++) {
				Eigen::VectorXd state = basis.col(k);
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

/* Fill the interaction term of the Hamiltonian */
void BH::fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double U) const {
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(hmatrix.nonZeros() + basis.cols());
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

/* Fill the chemical potential term of the Hamiltonian */
void BH::fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double mu) const {
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(hmatrix.nonZeros() + basis.cols());
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


    /// CONSTRUCTOR ///

/* Constructor for the Bose-Hubbard model */
BH::BH(const std::vector<std::vector<int>>& neighbours, int m, int n, double J, double U, double mu) : neighbours(neighbours), m(m), n(n), D(dimension(m,n)), J(J), U(U), mu(mu), H(D,D) {
    Eigen::MatrixXd basis = init_lexicographic(m, n);
    H.setZero();
    if (J != 0) {
        std::vector<int> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
        Eigen::VectorXd tags = calculate_tags(basis, primes);
        sort_basis(tags, basis);
        fill_hopping(basis, tags, neighbours, primes, H, J);
    }
    if (U != 0) {
        fill_interaction(basis, H, U);
    }
    if (mu != 0) {
        fill_chemical(basis, H, mu);
    }
}


    /// UTILITY FUNCTIONS ///

/* get the Hamiltonian matrix */
Eigen::SparseMatrix<double> BH::getHamiltonian() const {
	return H;
}

/* Display all the attributes of the Hamiltonian */
void BH::display_all() const{
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 
    displaySparseMatrix(H);
}