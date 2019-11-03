#include "Potential.h"
#include <iostream>

// Constructor initializes dist and forces vectors, and links the Atoms object
Potential::Potential(Atoms* a)
	: dist(a->getSize(), vector<double>(a->getSize(), 0)),
	forces(a->getSize(), vector<double>(3, 0))
{
	atoms = a;
}

// Destructor releases the memory of the internal vectors
Potential::~Potential() {
	vector<vector<double>>().swap(dist);
	vector<vector<double>>().swap(forces);
}

// Getter for the sumForceInteraction member
double Potential::getSumForcesInteraction() {
	return sumForceInteractions;
}

// Let the forces be recalculated
void Potential::reset() {
	forcesCalculated = false;
}

// Constructor for the Lennard-Jones potential initializes as a Potential
LJ::LJ(Atoms* a) :
	Potential(a) {}

// Function for returning the potential energy 
double LJ::getEnergy() {
	// Get the distances
	dist = atoms->getDistances();

	// Initialize the energy as zero
	double U = 0.0;
	// Run over all atom pairs and calculate the energy
	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			U += 4.0 * (pow(1.0 / dist[i][j], 12.0) - pow(1.0 / dist[i][j], 6.0));
		}
	}
	// Return the potential energy
	return U;
}

vector<vector<double>> LJ::getForces() {
	// Only recalculate the forces, if they need to be
	if (forcesCalculated) {
		return forces;
	}
	// Get the distances
	dist = atoms->getDistances();
	
	// Reset the sumForceInteractions
	sumForceInteractions = 0.0;

	// Run through all atom pairs
	vector<vector<double>> F(atoms->getSize(), vector<double>(3, 0));
	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			// scalar force divided by distance
			double pf = 48 / dist[i][j] * (pow(1 / dist[i][j], 14.0)
				- 0.5 * pow(1 / dist[i][j], 8.0));
			for (int k = 0; k < 3; k++)	{
				// Calculate the pbc distance per axis
				double diff = atoms->getPos(i)[k] - atoms->getPos(j)[k];
				double pbc_dist = diff - atoms->getCellLength()
					* round(diff / atoms->getCellLength());

				// Multiply the prefactor with the distance
				double F_jia = pf * pbc_dist;
				
				// Add the force to the vector of both affected atoms
				F[i][k] += F_jia;
				F[j][k] -= F_jia;
				// Add the force interaction
				sumForceInteractions += F_jia * pbc_dist;
			}
		}
	}
	// Save the forces to the internal memory
	forces = F;
	// The forces have now been calculated
	forcesCalculated = true;
	return F;
}

// Helper function for printing the distance matrix to the console
void LJ::printDistances() {
	for (vector<double> p : dist) {
		for (double d : p) {
			cout << d << ", ";
		}
		cout << "\n";
	}
}

// Helper function for printing the forces vector to the console
void LJ::printForces(vector<vector<double>> F) {
	for (vector<double> f : F) {
		for (double a : f) {
			cout << a << ", ";
		}
		cout << endl;
	}
}