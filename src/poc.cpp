#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<array>
#include<cmath>


#include<Eigen> 
// #include<xdiag>


#include"src/lattice_xy.h"


int main()
{

    Lattice1D_XY chain(3);
    chain.display_all();
    
    return 0; 
}
