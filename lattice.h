#ifndef LATTICE_H
#define LATTICE_H


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<array>

#include<Eigen>

class Lattice

/* Lattice class implements the lattice, be it 1D or 2D, ie the finite Hilbert Space of our
Quantum System. This class implements the model XY used to describe the system.
NB: Lattice is an abstract class, only its subcclasses Lattice1D or Lattice2D can be 
instantiated */  

{

    friend class Operator;   

    virtual void display_all() const = 0; //pure virtual function
    // function that displays all the informations on the lattice, and a 
    // visualization of the lattice if possible (use SFML ?) 

}; 


class Lattice1D : public Lattice
/* Implementation of a spin chain, using the XY model (quantum version)
A basis of the Hilbert Space of dimension D = 2^N is for instance 
composed of the tensor products states of the spin 1/2 states of each site.
All the operators are represented in that specific basis */
{
    private: 

    // Attributes

    // Initialization of the Pauli Matrices (but the declaration has to be done using a constructor), see lattie.cpp file
    Eigen::SparseMatrix<std::complex<double>> S_x; 
    Eigen::SparseMatrix<std::complex<double>> S_y;
    Eigen::SparseMatrix<std::complex<double>> S_z;  

    // Identity matrix
    Eigen::SparseMatrix<std::complex<double>> I; 


    unsigned int N; // number of sites, ie of spins, in the chain
    unsigned int D; // dimension of the Hilbert space of the chain

    
    Eigen::SparseMatrix<std::complex<double>> H; // sparse matrix representation of the Hamiltonian of the chain


    // Operator overloading

    
    //Methods

    Eigen::SparseMatrix<std::complex<double>> kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::SparseMatrix<std::complex<double>>& B) const; // function that computes the tensor product of two sparse matrices A and B

    std::array<unsigned int, 2> neighbours(unsigned int i) const; 

    public: 


    // Constructors and Destructor
    Lattice1D(unsigned int N_); // Constructs a Spin Chain of N sites
                                // The initial quantum state is by default the ground state of the system 

    Lattice1D(unsigned int N_, bool random_is_true); // Constructs a Spin Chain of N sites
                                                     // The initial quantum state is a random state
                                                     // superposition of some of the tensor product in the basis state 

    


    ~Lattice1D(); // Destructor


    // Display methods 

    void display_all() const; 
    void displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const; // function that displays a sparse matrix


};


class Lattice2D : public Lattice
/*Implementation of a 2D spin lattice, using the XY model (classical and quantum version) */
{

}; 


#endif // LATTICE_H