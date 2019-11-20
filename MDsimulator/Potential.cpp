#define _USE_MATH_DEFINES
#include "Potential.h"
#include <iostream>

// Constructor initializes dist and forces vectors, and links the Atoms object.
// If the radial cut-off is in use, it also calculates constants for this.
Potential::Potential(Atoms* a, double nDensity, double cutoff)
	: dist(a->getSize(), vector<double>(a->getSize(), 0)),
	forces(a->getSize(), vector<double>(3, 0))
{
	atoms = a;
	numberDensity = nDensity;
	r_c = cutoff;
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

// Constructor for the Lennard-Jones potential initializes as a Potential
LJ::LJ(Atoms* a, double nDensity, double cutoff) :
	Potential(a, nDensity, cutoff) 
{
	if (r_c != 0.0) {
		cutoffEnergy = calculateEnergy(r_c);
		diffU_r = -48 * (pow(1.0 / r_c, 13.0) - 0.5 * pow(1.0 / r_c, 7.0));
	}
}

// Function for returning the potential energy 
double LJ::getEnergy() {
	// Get the distances
	dist = atoms->getDistances();

	// Initialize the energy as zero
	double U = 0.0;
	// Run over all atom pairs and calculate the energy
	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			U += calculateEnergy(dist[i][j]);
		}
	}
	// Return the potential energy
	return U + calculateEnergyCorrection();
}

vector<vector<double>> LJ::getForces() {
	// Only recalculate the forces, if the atomic positions have changed
	if (!(atoms->hasChangedPositions())) {
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
			if (r_c != 0.0 && dist[i][j] > r_c) {
				continue;
			}
			// force prefactor
			double pf = 48 * ( pow(1.0 / dist[i][j], 14.0) 
				- 0.5 * pow(1.0 / dist[i][j], 8.0) );
			for (int k = 0; k < 3; k++)	{
				// Calculate the pbc distance per axis
				double diff = atoms->getPos(i)[k] - atoms->getPos(j)[k];
				double pbc_dist = diff - atoms->getCellLength()
					* round(diff / atoms->getCellLength());

				// Multiply the prefactor with the distance
				double F_jia = pf * pbc_dist;
				
				// If we are working with a cut-off, then add the correction
				if (r_c != 0.0) {
					F_jia += diffU_r * pbc_dist / dist[i][j];
				}

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
	return F;
}

double LJ::getPressureCorrection() {
	if (r_c == 0.0) {
		return 0;
	}
	return 32.0 / 9.0 * M_PI * numberDensity * numberDensity *
		(pow(1.0 / r_c, 9.0) - 1.5 * pow(1.0 / r_c, 3.0));
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

double LJ::calculateEnergy(double r) {
	if (r_c != 0.0 && r_c < r) {
		return 0;
	}
	double U_r = 4.0 * (pow(1.0 / r, 12.0) - pow(1.0 / r, 6.0));
	if (r_c == 0.0) {
		return U_r;
	}
	return U_r - cutoffEnergy - diffU_r * (r - r_c);
}

double LJ::calculateEnergyCorrection() {
	if (r_c == 0.0) {
		return 0;
	}
	return 8.0 / 9.0 * M_PI * atoms->getSize() * numberDensity *
		(pow(1.0 / r_c, 9.0) - 3.0 * pow(1.0 / r_c, 3.0));
}