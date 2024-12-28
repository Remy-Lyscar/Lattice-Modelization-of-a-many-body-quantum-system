#ifndef LATTICE_BH_H
#define LATTICE_BH_H


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

    int N; // number of sites in the chain
    int M; // number of bosons in the chain
    
    int D; // dimension of the Hilbert space of the chain
    double J; // Coupling constant of the XY model 
    double U; // Interaction parameter of the BH model 
    double mu; // tranverse magnetic field parameter of the XY model 

    Eigen::SparseMatrix<std::complex<double>> H; // sparse matrix representation of the Bose-Hubbard Hamiltonian of the chain


    // Private methods to build the basis vectors of the Hilbert space and then set the Bose-Hubbard Hamiltonian

    int factorial(int n) const;
    int dimension(int N_, int M_) const;
    int sum(const Eigen::VectorXd& state, int index1, int index2) const;
    bool next_lexicographic(Eigen::VectorXd& state, int m, int n) const;
    Eigen::MatrixXd basis_lexicographic(int m, int n) const;



    public: 

    // Constructors and Destructor
    Lattice1D_BH(int N_, int M_); // Default constuctor, constructs a chain of N sites and 
                                       // M bosons, with coupling parameter J = 1.0, interaction parameter U = 1.0 
                                       // and mu = 0.0


    Lattice1D_BH(int N_,int M_,  double J_, double U_);

    Lattice1D_BH(int N_, int M_, double J_, double U_, double mu_);     


    ~Lattice1D_BH(); // Destructor


    // Display methods 

    void display_all() const; 
    void displaySparseMatrix(const Eigen::SparseMatrix<std::complex<double>>& M) const; // function that displays a sparse matrix


}; 


#endif // LATTICE_BH