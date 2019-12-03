#include "Atoms.h"
#include <iostream>
#include "dataType.h"

// The constructor initializes the position and velocity vectors and the
// distance matrix to the right size
Atoms::Atoms(int natoms, double m)
	: pos(natoms, vector<double>(3, 0)),
	vel(natoms, vector<double>(3, 0)),
	distances(natoms, vector<double>(natoms, 0))
{
	// Initialize the number of atoms and the cell size
	nAtoms = natoms;
	cellLength = 0.0;
	mass = m;
}

// The destructor deletes the memory of the position and velocity vectors
Atoms::~Atoms() {
	vector<vector<double>>().swap(pos);
	vector<vector<double>>().swap(vel);
	vector<vector<double>>().swap(distances);
}

// The print() function print out all the positions, and an average position
void Atoms::print() {
	// Create average vector
	vector<double> R = { 0.0, 0.0, 0.0 };
	// Run through all positions
	for (vector<double> p : pos) {
		for (int i = 0; i < 3; i++) {
			R[i] += p[i];
			cout << p[i] << ", ";
		}
		cout << endl;
	}
	// Calculate and print the average
	cout << "Average: ";
	for (double d : R) {
		cout << d / nAtoms << ", ";
	}
	cout << endl;
}

// The printVel() function print out all the velocities, and an average velocity
void Atoms::printVel() {
	// Create average vector
	vector<double> av = { 0.0, 0.0, 0.0 };
	// Run through all velocities
	for (vector<double> v : vel) {
		for (int i = 0; i < 3; i++) {
			av[i] += v[i];
			cout << v[i] << ", ";
		}
		cout << "\n";
	}
	// Calculate and print the average
	cout << "Average: ";
	for (double d : av) {
		cout << d / nAtoms << ", ";
	}
	cout << endl;
}

// Helper function for printing the distance matrix to the console
void Atoms::printDistances() {
	if (positionsChanged) {
		distances = this->getDistances();
	}
	for (vector<double> p : distances) {
		for (double d : p) {
			cout << d << ", ";
		}
		cout << "\n";
	}
}

// The getDistances() function calculates all interatomic distances with
// periodic boundary conditions and returns these. The calculation is only
// performed, if the atomic positions have changed
vector<vector<double>> Atoms::getDistances() {
	// Don't recalculate, if the positions are the same
	if (!positionsChanged) {
		return distances;
	}

	// Run through each pair and calculate the distances
	for (int i = 0; i < nAtoms - 1; i++) {
		for (int j = i + 1; j < nAtoms; j++) {
			double r = 0.0;  // distance
			for (int k = 0; k < 3; k++) {
				double diff = pos[i][k] - pos[j][k];
				// Periodic Boundary Condition distance
				double pbc_dist = diff - cellLength	* round(diff / cellLength);
				// square the cartesian distance
				r += pbc_dist * pbc_dist;
			}
			r = pow(r, 0.5);  // euclidean space r = (x^2 + y^2 + z^2)^(1/2)
			// Add the distance to both sides of the matrix
			distances[i][j] = r;
			distances[j][i] = r;
		}
	}
	// The positions are now not changed
	positionsChanged = false;
	return distances;
}

// The center() function makes the average position of all the atoms (0, 0, 0)
void Atoms::center() {
	vector<double> R = { 0.0, 0.0, 0.0 };
	double M = 0;
	// find the center of mass (COM)
	for (vector<double> p : pos) {
		for (int i = 0; i < 3; i++) {
			R[i] += mass * p[i];
		}
		M += mass;
	}
	for (int i = 0; i < 3; i++) {
		R[i] /= M;
	}

	// Adjust all positions, so they lie around COM
	for (vector<double> &p : pos) {
		for (int i = 0; i < 3; i++)	{
			p[i] -= R[i];
		}
	}
}

// The centerVel() function makes the average velocity of all the atoms (0, 0, 0)
void Atoms::centerVel() {
	vector<double> av = { 0.0, 0.0, 0.0 };
	double M = 0.0;
	// Find the mass-weigthed center of velocity (COV)
	for (vector<double> v : vel) {
		for (int i = 0; i < 3; i++)	{
			av[i] += mass * v[i];
		}
		M += mass;
	}
	for (int i = 0; i < 3; i++) {
		av[i] /= M;
	}

	// Adjust all velocities, so they lie around COV
	for (vector<double> &v : vel) {
		for (int i = 0; i < 3; i++)	{
			v[i] -= av[i];
		}
	}
}

// Simple getter for the number of atoms
int Atoms::getSize() {
	return nAtoms;
}

// Simple getter for the cell size
double Atoms::getCellLength() {
	return cellLength;
}

// Simple getter for the mass
double Atoms::getMass() {
	return mass;
}

// Getter for the position vector of atom i
vector<double> Atoms::getPos(int i) {
	return pos[i];
}

// Getter for the velocity vector of atom i
vector<double> Atoms::getVel(int i) {
	return vel[i];
}

// The getEnergy() function returns the kinetic energy of all the atoms
double Atoms::getEnergy() {
	double K = 0.0;
	// Run over all cartesian coordinates of the atoms and add the velocity
	// squared to the kinetic energy
	for (vector<double> v : vel) {
		for (double vj : v) {
			K += vj * vj;
		}
	}
	// Multiply by the factor of a half
	return K / 2.0;
}

//
bool Atoms::hasChangedPositions() {
	return positionsChanged;
}

// Setter for the position vector of atom i
void Atoms::setPos(int i, vector<double> r) {
	positionsChanged = true;
	pos[i] = r;
}

// Setter for the velocity vector of atom i
void Atoms::setVel(int i, vector<double> r) {
	vel[i] = r;
}

// Setter for the cell size
void Atoms::setCellLength(double length) {
	// Ensure that the length is a positive number
	if (length > 0.0) {
		cellLength = length;
	}
}