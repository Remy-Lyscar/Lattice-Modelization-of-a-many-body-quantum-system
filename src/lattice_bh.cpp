#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<Eigen/Dense> 
#include<Eigen/SparseCore>
#include<complex>
#include<array>

#include"lattice_bh.h"



/*-----     Constructors and Destructor     -----*/

Lattice1D_BH::Lattice1D_BH(int N_, int M_): N(N_), M(M_), D(dimension(N_, M_)), J(1.0), U(1.0), mu(0.0), H(D,D), B(D,D)
{
    // Initialization of the Hamiltonian of the Bose-Hubbard model
    //computeHamiltonianBH(); 

    B = basis_lexicographic(N, M); // Basis vectors of the Hilbert space of the chain

    display_matrix(B); // Display the basis vectors of the Hilbert space of the chain
    // initialization of the initial state of the Bose-Hubbard chain, ie the ground state of the system

}



Lattice1D_BH::~Lattice1D_BH()
{
    std::cout << "The lattice implementing Bose-Hubbard model has been destroyed!" << std::endl; 
}


/*-----    Private: Utility methods     -----*/

int Lattice1D_BH::factorial(int N_) const
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

int Lattice1D_BH::dimension(int N_, int M_) const
{
    return (factorial(N_ + M_ - 1))/(factorial(N_)*factorial(M_ - 1)); 
}



/* calculate the sum of the elements of a vector between 2 index */
int Lattice1D_BH::sum(const Eigen::VectorXd& state, int index1, int index2) const { 
	int s = 0;
	for (int i = index1; i <= index2; i++) {
		s += state[i];
	}
	return s;
}


// INITIALIZE THE HILBERT SPACE BASIS

/* calculate the next state of the Hilbert space in lexicographic order */
bool Lattice1D_BH::next_lexicographic(Eigen::VectorXd& state, int m, int n) const {
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
Eigen::MatrixXd Lattice1D_BH::basis_lexicographic(int m, int n) const {
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

void Lattice1D_BH::display_all() const
{
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 

    displaySparseMatrix(H);
}


void Lattice1D_BH::displaySparseMatrix(const Eigen::SparseMatrix<double>& M) const
/* Display a sparse matrix in a standard form in the terminal, as a full matrix
The SparseMatrix given in argument is passed by reference to avoid the copy of the matrix 
in the local context of this method */
{
    std::cout << "Displaying the sparse matrix at memory adress: " << &M << "\n \n" << std::endl; 

    for(int i = 0; i < M.rows(); i++)
    {
        for(int j = 0; j < M.cols(); j++)
        {
            std::complex<double> c = M.coeff(i, j);
            if (std::norm(c) < epsilon)
            {
                std::cout << "(0,0)" << " ";  // We consider that the coefficient is indeed zero 
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


void Lattice1D_BH::display_matrix(const Eigen::MatrixXd& M) const
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



