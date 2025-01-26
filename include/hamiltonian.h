#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include<vector>
#include<complex>

#include<Eigen/Dense>
#include<Eigen/SparseCore> 




/**
 * @brief Base class for Hamiltonians.
 * 
 * This class provides a common interface for different types of Hamiltonians.
 */
class Hamiltonian{
public :

// NEIGHBOURS

    /**
     * @brief Generate the list of neighbours for a 1D chain.
     * 
     * @param m Number of sites in the chain.
     * @param closed By default, closed = true for periodic boundary conditions, closed = false for open boundary conditions.
     * @return std::vector<std::vector<int>> The list of neighbours for each site of the chain.
     */
    std::vector<std::vector<int>> chain_neighbours(int m, bool closed = true);

    /**
     * @brief Generate the list of neighbours for a 2D square lattice.
     * 
     * @param m Number of sites in the square.
     * @param closed By default, closed = true for periodic boundary conditions, closed = false for open boundary conditions.
     * @return std::vector<std::vector<int>> The list of neighbours for each site of the square.
     */
     std::vector<std::vector<int>> square_neighbours(int m, bool closed = true);

// DISPLAY METHODS

    /**
     * @brief Display all the attributes of the Hamiltonian.
     */
    virtual void display_all() const = 0; 
    
    /**
     * @brief Display a sparse matrix in a standard form in the terminal, as a full matrix.
     * 
     * @tparam T The type of the elements in the sparse matrix.
     * @param M The sparse matrix to display.
     */
    template <typename T> void displaySparseMatrix(const Eigen::SparseMatrix<T>& M) const; 

    /**
     * @brief Display a dense matrix in the terminal.
     * 
     * @param M The dense matrix to display.
     */
    void displayMatrix(const Eigen::MatrixXd& M) const; 

protected :
    double epsilon = 1e-10; // threshold beneath which a number is considered as zero 

}; 




/**
 * @brief Class representing the XY Hamiltonian.
 * 
 * This class implements the Hamiltonian for the XY model, which is a quantum spin chain model.
 */
class XY : public Hamiltonian{
private : 

// PARAMETERS OF THE XY MODEL

    int N; // number of sites, ie of spins, in the chain
    int D; // dimension of the Hilbert space of the chain

    double J; // Coupling constant of the XY model 
    double mu; // Tranverse magnetic field parameter of the XY model 

    Eigen::SparseMatrix<std::complex<double>> H; // Sparse matrix representation of the Hamiltonian 

// INITIALIZATION OF THE PAULI MATRICES AND THE IDENTITY MATRIX

    Eigen::SparseMatrix<std::complex<double>> S_x; 
    Eigen::SparseMatrix<std::complex<double>> S_y;
    Eigen::SparseMatrix<std::complex<double>> S_z;  
    Eigen::SparseMatrix<std::complex<double>> I;

//METHODS

    /* function that computes the Kronecker product of two sparse matrices */
    Eigen::SparseMatrix<std::complex<double>> kroneckerProductSparse(const Eigen::SparseMatrix<std::complex<double>>& A,
                                                                     const Eigen::SparseMatrix<std::complex<double>>& B) const; 

    /* function that computes the Kronecker product of a Pauli matrix S at site i with the identity matrix at the other sites */
    Eigen::SparseMatrix<std::complex<double>> kroneckerPauli(const Eigen::SparseMatrix<std::complex<double>>& S, int i) const; 
    void computeHamiltonianXY(); 
    
public: 

// CONSTRUCTORS

    /**
     * @brief Constructor for the XY Hamiltonian with by default coupling parameter J = 1.0 and initial quantum state the ground state of the system
     * 
     * @param N Number of sites in the chain.
     */
    XY(int N_); 

    /**
     * @brief Constructor for the XY Hamiltonian with coupling parameter J.
     * 
     * @param N Number of sites in the chain.
     * @param J Coupling parameter of the XY model.
     */
    XY(int N_, double J_); 

    /**
     * @brief Constructor for the XY Hamiltonian with coupling parameter J and transverse magnetic field parameter mu.
     * 
     * @param N Number of sites in the chain.
     * @param J Coupling parameter of the XY model.
     * @param mu Transverse magnetic field parameter of the XY model.
     */
    XY(int N_, double J_, double mu_);

    /**
     * @brief Constructor for the XY Hamiltonian with a random initial state.
     * 
     * @param N Number of sites in the chain.
     * @param random_is_true If true, the initial quantum state is a random state superposition of some of the tensor product in the basis state.
     */
    XY(int N_, bool random_is_true); 
                

// UTILITY FUNCTIONS

    /**
    * @brief Return the Hamiltonian sparse matrix.
    * 
    * @return Eigen::SparseMatrix<std::complex<double>> The Hamiltonian sparse matrix.
    */
    Eigen::SparseMatrix<std::complex<double>> getHamiltonian() const; 

// DISPLAY METHODS

    /**
    * @brief Display all the attributes of the XY model.
    */
    void display_all() const; 
}; 




/**
 * @brief Class representing the Bose-Hubbard Hamiltonian.
 * 
 * This class implements the Hamiltonian for the Bose-Hubbard model, which describes interacting bosons on a lattice.
 */
class BH: public Hamiltonian{
private:

// PARAMETERS OF THE BH MODEL

    std::vector<std::vector<int>> neighbours; // Vector that contains the neighbours of each site of the chain
    
    int m; // Number of sites in the chain
    int n; // Number of bosons in the chain
    int D; // Dimension of the Hilbert space of the system

    double J; // Hopping parameter of the BH model
    double U; // Interaction parameter of the BH model 
    double mu; // Chemical potential of the BH model

    Eigen::SparseMatrix<double> H; // Sparse matrix representation of the Hamiltonian

// DIMENSION OF THE HILBERT SPACE

    /* calculate the dimension of the Hilbert space for n bosons on m sites */
    int binomial(int n, int k) const; // Binomial coefficient
    int dimension(int m, int n) const; // Dimension of the Hilbert space

// ELEMENTARY FUNCTIONS

    /* calculate the sum of the elements of a vector between 2 index */
    int sum(const Eigen::VectorXd& state, int index1, int index2) const; 

// INITIALIZE THE HILBERT SPACE BASIS

    /* calculate the next state of the Hilbert space in lexicographic order */
    bool next_lexicographic(Eigen::VectorXd& state, int m, int n) const; 

    /* creates a matrix that has the vectors of the Hilbert space basis in columns */
    Eigen::MatrixXd init_lexicographic(int m, int n) const; 

// SORT THE HILBERT SPACE BASIS TO FACILITATE CALCULUS

    /* calculate the unique tag of the kth column of the matrix */
    double calculate_tag(const Eigen::MatrixXd& basis, const std::vector<int>& primes, int k) const;

    /* calculate and store the tags of each state of the Hilbert space basis */
    Eigen::VectorXd calculate_tags(const Eigen::MatrixXd& basis, const std::vector<int>& primes) const;

    /* sort the states of the Hilbert space by ascending order compared by their tags */
    void sort_basis(Eigen::VectorXd& tags, Eigen::MatrixXd& basis) const; 

    /* gives the index of the wanted tag x by the Newton method */
    int search_tag(const Eigen::VectorXd& tags, double x) const;

// FILL THE HAMILTONIAN OF THE SYSTEM

    /* fill the hopping term of the Hamiltonian */
    void fill_hopping(const Eigen::MatrixXd& basis, const Eigen::VectorXd& tags, const std::vector<std::vector<int>>& neighbours, const std::vector<int>& primes, Eigen::SparseMatrix<double>& hmatrix, double J) const;

    /* fill the interaction term of the Hamiltonian */
    void fill_interaction(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double U) const; 

    /* fill the chemical potential term of the Hamiltonian */
    void fill_chemical(const Eigen::MatrixXd& basis, Eigen::SparseMatrix<double>& hmatrix, double mu) const; 

public:

// CONSTRUCTOR

    /**
    * @brief Constructor for the Bose-Hubbard Hamiltonian.
    * 
    * @param neighbours Vector that contains the neighbours of each site of the lattice.
    * @param m Number of sites in the lattice.
    * @param n Number of bosons in the lattice.
    * @param J Hopping parameter of the BH model.
    * @param U Interaction parameter of the BH model.
    * @param mu Chemical potential of the BH model.
    */
    BH(const std::vector<std::vector<int>>& neighbours, int m_, int n_, double J_, double U_, double mu_);

// UTILITY FUNCTIONS

    /**
    * @brief Return the Hamiltonian sparse matrix.
    * 
    * @return Eigen::SparseMatrix<double> The Hamiltonian sparse matrix.
    */
    Eigen::SparseMatrix<double> getHamiltonian() const; 

// DISPLAY METHODS

    /**
    * @brief Display all the attributes of the BH model.
    */
    void display_all() const; 
};

#endif // HAMILTONIAN_H
 