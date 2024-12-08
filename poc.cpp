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


#include"lattice.h"


int main()
{

    Lattice1D chain(3);
    chain.display_all();
    
    return 0; 
}
