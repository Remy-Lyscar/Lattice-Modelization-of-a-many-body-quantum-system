#ifndef LATTICE
#define LATTICE


#include<chrono>
#include<random>
#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<array>

#include<Eigen/Dense>
#include<Eigen/SparseCore> 

class Lattice
{
    public : 

    virtual ~Lattice() = default;
    virtual void display_all() const = 0; 
    
    template <typename T> void displaySparseMatrix(const Eigen::SparseMatrix<T>& M) const; // function that displays a sparse matrix in a standard form in the terminal, as a full matrix
    void display_matrix(const Eigen::MatrixXd& M) const; // function that displays a dense matrix

    protected :

    double epsilon = 1e-10; // threshold beneath which a number is considered as zero

    std::array<int, 2> neighbours(int i) const; 

}; 

class XY : public Lattice
{

    private : 
    
    // Attributes

    int N; // number of sites, ie of spins, in the chain
    
    int D; // dimension of the Hilbert space of the chain
    double J; // Coupling constant of the XY model 
    double mu; // tranverse magnetic field parameter of the XY model 


    // Initialization of the Pauli Matrices 
    Eigen::SparseMatrix<std::complex<double>> S_x; 
    Eigen::SparseMatrix<std::complex<double>> S_y;
    Eigen::SparseMatrix<std::complex<double>> S_z;  

    // Identity matrix
    Eigen::SparseMatrix<std::complex<double>> I; 


    // Operator overloading

    
    //Methods

    Eigen::SparseMatrix<std::complex<double>> kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A,
                                                                     const Eigen::SparseMatrix<std::complex<double>>& B) const; // function that computes the tensor product of two sparse matrices A and B
    
    Eigen::SparseMatrix<std::complex<double>> kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, int i) const; 
    void computeHamiltonianXY(); // function that computes the Hamiltonian of the XY model
    
    std::array<int, 2> neighbours(int i) const; 

    public: 

    // Attributes
    Eigen::SparseMatrix<std::complex<double>> H; // sparse matrix representation of the Hamiltonian of the chain

    // Constructors and Destructor
    XY(int N_); // Constructs a Spin Chain of N sites, with coupling parameter J = 1.0
                                // The initial quantum state is by default the ground state of the system 


    XY(int N_, double J_); 

    XY(int N_, double J_, double mu_);

    XY(int N_, bool random_is_true); // Constructs a Spin Chain of N sites
                                    // The initial quantum state is a random state
                                    // superposition of some of the tensor product in the basis state 

    


    ~XY(); // Destructor


    // Display methods 

    void display_all() const; 
}; 




class BH: public Lattice
{
    private: 

    // Attributes

    int N; // number of sites in the chain
    int M; // number of bosons in the chain
    
    int D; // dimension of the Hilbert space of the chain
    double J; // Coupling constant of the XY model 
    double U; // Interaction parameter of the BH model 
    double mu; // tranverse magnetic field parameter of the XY model 

    Eigen::SparseMatrix<double> H; // sparse matrix representation of the Bose-Hubbard Hamiltonian of the chain
    Eigen::MatrixXd B; // matrix that contains the basis vectors of the Hilbert space of the chain

    // Private methods to build the basis vectors of the Hilbert space and then set the Bose-Hubbard Hamiltonian

    int factorial(int n) const;
    int dimension(int N_, int M_) const;
    int sum(const Eigen::VectorXd& state, int index1, int index2) const;
    bool next_lexicographic(Eigen::VectorXd& state, int m, int n) const;
    Eigen::MatrixXd basis_lexicographic(int m, int n) const;

    // MatrixXd is an alias for: Matrix<double, Dynamic, Dynamic>



    public: 

    // Constructors and Destructor
    BH(int N_, int M_); // Default constuctor, constructs a chain of N sites and 
                                       // M bosons, with coupling parameter J = 1.0, interaction parameter U = 1.0 
                                       // and mu = 0.0


    BH(int N_,int M_,  double J_, double U_);

    BH(int N_, int M_, double J_, double U_, double mu_);     


    ~BH(); // Destructor


    // Display methods 

    void display_all() const; 
}; 


#endif  // LATTICE