#ifndef LATTICE_H
#define LATTICE_H


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>

#include<Eigen>

class Lattice

/* Lattice class implements the lattice, be it 1D or 2D, ie the finite Hilbert Space of our
Quantum System. This class implements the model XY used to describe the system.
NB: Lattice is an abstract class, only its subcclasses Lattice1D or Lattice2D can be 
instantiated */  

{

    friend class Operator;   

    virtual unsigned int memory() const = 0; //pure virtual function
    // function that returns the memory used by the lattice, in octets 

}; 


class Lattice1D : public Lattice
/* Implementation of a spin chain, using the XY model (classical and quantum version)*/
{
    private: 

    // Initialization of the Pauli Matrices (but the declaration has to be done using a constructor), see lattie.cpp file
    Eigen::SparseMatrix<std::complex<double>> S_x; 
    Eigen::SparseMatrix<std::complex<double>> S_y;
    Eigen::SparseMatrix<std::complex<double>> S_z;  


    unsigned int N; // number of sites, ie of spins, in the chain
    unsigned int D; // dimension of the Hilbert space of the chain

    
    Eigen::SparseMatrix<std::complex<double>> H; // sparse matrix representation of the Hamiltonian of the chain


    Eigen::SparseMatrix<std::complex<double>> tensor_product(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::SparseMatrix<std::complex<double>>& B) const;
    // function that computes the tensor product of two sparse matrices A and B

    public: 

    Lattice1D(unsigned int N_); // Constructs a Spin Chain of N sites
    ~Lattice1D(); // Destructor



};


class Lattice2D : public Lattice
/*Implementation of a 2D spin lattice, using the XY model (classical and quantum version) */
{

}; 


#endif // LATTICE_H