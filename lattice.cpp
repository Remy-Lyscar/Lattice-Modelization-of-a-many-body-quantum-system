#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<Eigen> 
#include<complex>
#include<array>

#include"lattice.h"
// #include"operator.h"



/*-----     Constructors and Destructor     -----*/

Lattice1D::Lattice1D(unsigned int N_): N(N_), D(1<<N_), J(1.0), mu(0.0), S_x(2,2), S_y(2,2), S_z(2,2), I(2,2), H(D,D)
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

    // initialization of the initial state of the spin 1/2 chain, ie the ground state of the system


}




Lattice1D::~Lattice1D()
{
    std::cout << "Lattice1D Destructor called" << std::endl; 
}



/*-----     Private methods     -----*/


std::array<unsigned int, 2> Lattice1D::neighbours(unsigned int i) const
/* Returns the array of the indexes of the neihbours of a site i in the chain*/
{  
    std::array<unsigned int, 2> nei; 

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

Eigen::SparseMatrix<std::complex<double>> Lattice1D::kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::SparseMatrix<std::complex<double>>& B) const
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


Eigen::SparseMatrix<std::complex<double>> Lattice1D::kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, unsigned int i) const
/* Method that computes the Kronecker Product of a Pauli matrix S at site i with the identity matrix at the other sites */
{
    Eigen::SparseMatrix<std::complex<double>> result = (i==0) ? S: I;  

    for(unsigned int j = 1; j < N; j++)
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


void Lattice1D::computeHamiltonianXY() 
/* This method computes the full hamiltonian of the XY quantum system without keeping 
in memory the individual spin operators, that are useless for what is to come 
That method is not const since it's designed to modify class attribute H */
{

    // Initialize H with the first term. Operator * is overloaded for SparseMatrix multiplication
    H = kroneckerPauli(S_x,0)*kroneckerPauli(S_x, 1) + kroneckerPauli(S_y,0)*kroneckerPauli(S_y, 1);

    for(unsigned int i = 1; i < N-1; i++)
    {
        H += kroneckerPauli(S_x,i)*kroneckerPauli(S_x, i+1) + kroneckerPauli(S_y,i)*kroneckerPauli(S_y, i+1);
    }

}




/*-----     Public methods     -----*/

void Lattice1D::display_all() const
{
    std::cout << "Soon, you will see all the informations about the lattice here!" << std::endl; 

    displaySparseMatrix(H);
}


void Lattice1D::displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const
/* Display a sparse matrix in a standard form in the terminal, as a full matrix
The SparseMatrix given in argument is passed by referene to avoid the copy of the matrix 
in the local context of this method */
{
    std::cout << "Displaying the sparse matrix at memory adress: " << &M << "\n \n" << std::endl; 

    for(int i = 0; i < M.rows(); i++)
    {
        for(int j = 0; j < M.cols(); j++)
        {
            std::cout << M.coeff(i, j) << " "; 
        }
        std::cout << std::endl; 
    }

    std::cout << "\n \n" << std::endl;
}
