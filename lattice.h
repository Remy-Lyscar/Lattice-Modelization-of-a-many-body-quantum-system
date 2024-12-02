#ifndef LATTICE_H
#define LATTICE_H


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>

#include<Eigen>

class Latice

/* Lattice class implements the lattice, ie the finite Hilbert Space of our
Quantum System. This class implements the Bose-Hubbard Hamiltonian of the system */

{

    friend class Operator;  
    friend class Dynamics;  
    friend class Quantum_States; 

    private: 

    unsigned short int d; //dimension of the lattice, usually 1 or 2
    unsigned int m, n; // Number of sites, number of bosons, 
    unsigned int D; // dimension of the Hilbert Space



    public: 


    





 

    


}; 


#endif // LATTICE_H