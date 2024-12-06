#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<Eigen> 
#include<complex>


#include"lattice.h"
// #include"operator.h"


Lattice1D::Lattice1D(unsigned int N_): N(N_), D(1<<N_)
{

    using namespace std::complex_literals; 

    S_x.insert(0, 0) = 0.; 
    S_x.insert(1,0) = 1.; 
    S_x.insert(0,1) = 1.; 
    S_x.insert(1,1) = 0.;

    S_y.insert(0, 0) = 0.; 
    S_y.insert(1,0) = 0. + 1i; 
    S_y.insert(0,1) = 0. - 1i; 
    S_y.insert(1,1) = 0.; 

    S_z.insert(0, 0) = 1.; 
    S_z.insert(1,0) = 0.; 
    S_z.insert(0,1) = 0.; 
    S_z.insert(1,1) = -1.;



}




Lattice1D::~Lattice1D()
{
    // std::cout << "Lattice1D Destructor called" << std::endl; 
}





