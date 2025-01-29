#include <vector>
#include <stdexcept>
#include <cmath>

#include "neighbours.h"


///// IMPLEMENTATION OF THE NEIGHBOURS CLASS METHODS /////


/* Constructor for the Neighbours class */
Neighbours::Neighbours(int m) :  m(m), neighbours(m) {}

/*Destructor for the Neighbours class */
Neighbours::~Neighbours() {}

/* generate the list of neighbours for a 1D chain */
void Neighbours::chain_neighbours(bool closed) { // by default closed = true for periodic boundary conditions, closed = false for open boundary conditions
	for (int i = 0; i < m; ++i) {
		if (i > 0) {
			neighbours[i].push_back(i - 1); // Left neighbour
		}
		if (i < m - 1) {
			neighbours[i].push_back(i + 1); // Right neighbour
		}
	}
	if (closed) { // Periodic boundary conditions
		neighbours[0].push_back(m - 1); 
		neighbours[m - 1].push_back(0); 
	}
}

/* generate the list of neighbours for a 2D square lattice */
void Neighbours::square_neighbours(bool closed) { // by default closed = true for periodic boundary conditions, closed = false for open boundary conditions
    int side = static_cast<int>(std::sqrt(m));
    if (side * side != m) {
        throw std::invalid_argument("The number of sites (m) must be a perfect square.");
    }
    std::vector<std::vector<int>> neighbours(m);
    for (int i = 0; i < side; ++i) {
        for (int j = 0; j < side; ++j) {
            int index = i * side + j;
            if (j > 0) {
                neighbours[index].push_back(index - 1); // Left neighbour
            } else if (closed) {
                neighbours[index].push_back(index + side - 1);
            }
            if (j < side - 1) {
                neighbours[index].push_back(index + 1); // Right neighbour
            } else if (closed) {
                neighbours[index].push_back(index - side + 1);
            }
            if (i > 0) {
                neighbours[index].push_back(index - side); // Top neighbour
            } else if (closed) {
                neighbours[index].push_back(index + (side - 1) * side);
            }
            if (i < side - 1) {
                neighbours[index].push_back(index + side); // Bottom neighbour
            } else if (closed) {
                neighbours[index].push_back(index - (side - 1) * side);
            }
        }
    }
}


/* Get the list of neighbours */
std::vector<std::vector<int>> Neighbours::getNeighbours() const {
    return neighbours;
}