#ifndef LATTICE_XY_H
#define LATTICE_XY_H


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<array>

#include<Eigen>



class Lattice1D_XY
{

    private : 
    
    // Attributes

    unsigned int N; // number of sites, ie of spins, in the chain
    
    int D; // dimension of the Hilbert space of the chain
    double J; // Coupling constant of the XY model (we impose that it is a positive real parameter)
    double mu; // tranverse magnetic field parameter of the XY model (we impose that it is a positive real parameter)

    // Initialization of the Pauli Matrices (but the declaration has to be done using a constructor), see lattie.cpp file
    Eigen::SparseMatrix<std::complex<double>> S_x; 
    Eigen::SparseMatrix<std::complex<double>> S_y;
    Eigen::SparseMatrix<std::complex<double>> S_z;  

    // Identity matrix
    Eigen::SparseMatrix<std::complex<double>> I; 

    
    Eigen::SparseMatrix<std::complex<double>> H; // sparse matrix representation of the Hamiltonian of the chain


    // Operator overloading

    
    //Methods

    Eigen::SparseMatrix<std::complex<double>> kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::SparseMatrix<std::complex<double>>& B) const; // function that computes the tensor product of two sparse matrices A and B
    Eigen::SparseMatrix<std::complex<double>> kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, unsigned int i) const; 
    void computeHamiltonianXY(); // function that computes the Hamiltonian of the XY model
    
    std::array<unsigned int, 2> neighbours(unsigned int i) const; 

    public: 


    // Constructors and Destructor
    Lattice1D_XY(unsigned int N_); // Constructs a Spin Chain of N sites, with coupling parameter J = 1.0
                                // The initial quantum state is by default the ground state of the system 


    Lattice1D_XY(unsigned int N_, double J_); 

    Lattice1D_XY(unsigned int N_, bool random_is_true); // Constructs a Spin Chain of N sites
                                                     // The initial quantum state is a random state
                                                     // superposition of some of the tensor product in the basis state 

    


    ~Lattice1D_XY(); // Destructor


    // Display methods 

    void display_all() const; 
    void displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const; // function that displays a sparse matrix in a standard form in the terminal, as a full matrix

}; 



#endif  // LATTICE_XY_H