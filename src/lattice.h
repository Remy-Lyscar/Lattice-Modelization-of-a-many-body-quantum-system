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

    friend class Lanczos;   

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

    virtual void display_all() const = 0; //pure virtual function
    

};



class Lattice2D : public Lattice
/*Implementation of a 2D spin lattice, using the XY model (classical and quantum version) */
{

    virtual void display_all() const = 0; //pure virtual function

}; 



#endif // LATTICE_H