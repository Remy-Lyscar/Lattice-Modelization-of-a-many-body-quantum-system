#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/SparseCore>
#include<complex>
#include<array>

#include "arpack++/arlsmat.h" // ARPACK++ classes for matix operations
#include "arpack++/arlssym.h" // ARPACK++ classes for symmetric standard eigenvalue problems

#include "lattice.h"


/*-----     Implementation of Lattice class methods     -----*/


template <typename T> void Lattice::displaySparseMatrix(const Eigen::SparseMatrix<T>& M) const
/* Display a sparse matrix in a standard form in the terminal, as a full matrix
The SparseMatrix given in argument is passed by reference to avoid the copy of the matrix 
in the local context of this method */
{
    std::cout << "Displaying the sparse matrix at memory adress: " << &M << "\n \n" << std::endl; 

    for(int i = 0; i < M.rows(); i++)
    {
        for(int j = 0; j < M.cols(); j++)
        {
            T c = M.coeff(i, j);
            if (std::norm(c) < epsilon)
            {
                std::cout << "(0,0)" << " "; 
            }
            else
            {
                std::cout  << "\033[34m" << c << " "; // Here we print the non zero elements in blue 
                                                      // so theyr are easy to see in the terminal
                std::cout << "\033[0m"; // Reset the color to the default one
            } 
        }
        std::cout << std::endl; 
    }

    std::cout << "\n \n" << std::endl;
}

// Explicit template instantiation for the types you need
template void Lattice::displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const;
template void Lattice::displaySparseMatrix(const Eigen::SparseMatrix<double>& M) const;


void Lattice::display_matrix(const Eigen::MatrixXd& M) const
{
    std::cout << "Displaying the matrix at memory adress: " << &M << "\n \n" << std::endl; 

    for(int i = 0; i < M.rows(); i++)
    {
        for(int j = 0; j < M.cols(); j++)
        {
            std::cout << M(i, j) << " "; 
        }
        std::cout << std::endl; 
    }

    std::cout << "\n \n" << std::endl;
}


std::array<int, 2> XY::neighbours(int i) const
/* Returns the array of the indexes of the neighbours of a site i in the chain*/
{  
    std::array<int, 2> nei; 

    //Periodic boundary conditions 
    if(i == 0)
    {
         nei[0] = N-1; 
         nei[1] = 1; 
    }
    else if(i == N-1)
    {
         nei[0] = N-2; 
         nei[1] = 0; 
    }
    else
    {
         nei[0] = i-1; 
         nei[1] = i+1; 
    }

    return nei; 
}



/*-----     Implementation of XY class methods     -----*/


/*-----     Constructors and Destructor     -----*/

XY::XY(int N_): N(N_), D(1<<N_), J(1.0), mu(0.0), S_x(2,2), S_y(2,2), S_z(2,2), I(2,2), H(D,D)
{

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
    S_y = 0.5*S_y;

    S_z.insert(0, 0) = 1.; 
    S_z.insert(1,0) = 0.; 
    S_z.insert(0,1) = 0.; 
    S_z.insert(1,1) = -1.;
    S_z = 0.5*S_z;


    // Initialization of the Hamiltonian of the XY model
    computeHamiltonianXY(); 

    displaySparseMatrix(H); // Display the Hamiltonian of the XY model

}




XY::~XY()
{
    std::cout << "XY Destructor called" << std::endl; 
}



/*-----     Private methods     -----*/


Eigen::SparseMatrix<std::complex<double>> XY::kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A,
                                                                               const Eigen::SparseMatrix<std::complex<double>>& B) const
/* Method that computes the Kronecker Product of two Eigen::SparseMatrix objects
and stores the result in a Eigen::SparseMAtrix object*/
{
    typedef Eigen::SparseMatrix<std::complex<double>> SparseMatrix;
    typedef Eigen::Triplet<std::complex<double>> Triplet; // The Triplet class is a utility class used to build SparseMatrix easily
    std::vector<Triplet> tripletList;


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


Eigen::SparseMatrix<std::complex<double>> XY::kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, int i) const
/* Method that computes the Kronecker Product of a Pauli matrix S at site i with the identity matrix at the other sites */
{
    Eigen::SparseMatrix<std::complex<double>> result = (i==0) ? S: I;  

    for(int j = 1; j < N; j++)
    {
        if(j == i)
        {
            result = kroneckerProductSparse(result, S); 
        }
        else
        {
            result = kroneckerProductSparse(result, I); 
        }
    }

    return result; 
}


void XY::computeHamiltonianXY() 
/* This method computes the full hamiltonian of the XY quantum system without keeping 
in memory the individual spin operators, that are useless for what is to come 
That method is not const since it's designed to modify class attribute H */
{

    // Initialize H with the first term. Operator * is overloaded for SparseMatrix multiplication
    H = kroneckerPauli(S_x,0)*kroneckerPauli(S_x, 1) + kroneckerPauli(S_y,0)*kroneckerPauli(S_y, 1);

    for(int i = 1; i < N-1; i++)
    {
        H += kroneckerPauli(S_x,i)*kroneckerPauli(S_x, i+1) + kroneckerPauli(S_y,i)*kroneckerPauli(S_y, i+1);
    }

}




/*-----     Public methods     -----*/

void XY::display_all() const
{
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 

    displaySparseMatrix(H);
}




/*-----     Implementation of Bose-Hubabrd class methods -----*/

/*-----     Constructors and Destructor     -----*/

BH::BH(int N_, int M_): N(N_), M(M_), D(dimension(N_, M_)), J(1.0), U(1.0), mu(0.0), H(D,D), B(D,D)
{
    // Initialization of the Hamiltonian of the Bose-Hubbard model
    //computeHamiltonianBH(); 

    B = basis_lexicographic(N, M); // Basis vectors of the Hilbert space of the chain

    display_matrix(B); // Display the basis vectors of the Hilbert space of the chain

    
    // initialization of the initial state of the Bose-Hubbard chain, ie the ground state of the system

}



BH::~BH()
{
    std::cout << "The lattice implementing Bose-Hubbard model has been destroyed!" << std::endl; 
}


/*-----    Private: Utility methods     -----*/

int BH::factorial(int N_) const
{
    if(N_ == 0)
    {
        return 1; 
    }
    else
    {
        return N_*factorial(N_-1); 
    }
}

int BH::dimension(int N_, int M_) const
{
    return (factorial(N_ + M_ - 1))/(factorial(N_)*factorial(M_ - 1)); 
}



/* calculate the sum of the elements of a vector between 2 index */
int BH::sum(const Eigen::VectorXd& state, int index1, int index2) const { 
	int s = 0;
	for (int i = index1; i <= index2; i++) {
		s += state[i];
	}
	return s;
}


// INITIALIZE THE HILBERT SPACE BASIS

/* calculate the next state of the Hilbert space in lexicographic order */
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

/* creates a matrix that has the vectors of the Hilbert space basis in columns
Rq: we actually construct the transpose */
Eigen::MatrixXd BH::basis_lexicographic(int m, int n) const {
	Eigen::MatrixXd basis(m, 1);
	basis.col(0).setZero();
	basis(0, 0) = n;  // |N, 0, 0, ..., 0> is the higer state of the Hilbert space in the lexicographic order

	Eigen::VectorXd state = basis.col(0);

	while (next_lexicographic(state, m, n)) {
		basis.conservativeResize(Eigen::NoChange, basis.cols() + 1); // Resize the matrix while conserving its current data 
		basis.col(basis.cols() - 1) = state;
	}
	return basis;
}


/*-----     Public methods     -----*/

void BH::display_all() const
{
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 

    displaySparseMatrix(H);
}


