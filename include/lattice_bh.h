#ifndef LATTICE_XY_H
#define LATTICE_XY_H


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<array>

#include<Eigen/Dense>
#include<Eigen/SparseCore> 


class Lattice1D_BH
{
    private: 

    // Attributes

    unsigned int N; // number of sites, ie of spins, in the chain
    
    int D; // dimension of the Hilbert space of the chain
    double J; // Coupling constant of the XY model (we impose that it is a positive real parameter)
    double mu; // tranverse magnetic field parameter of the XY model (we impose that it is a positive real parameter)

    public: 


    // Constructors and Destructor
    Lattice1D_BH(unsigned int N_); // Constructs a Spin Chain of N sites, with coupling parameter J = 1.0
                                // The initial quantum state is by default the ground state of the system 


    Lattice1D_BH(unsigned int N_, double J_);

    Lattice1D_BH(unsigned int N_, double J_, double mu_);     


    ~Lattice1D_BH(); // Destructor


    // Display methods 

    void display_all() const; 
    void displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const; // function that displays a sparse matrix


}; 


#endif // LATTICE_BH_H