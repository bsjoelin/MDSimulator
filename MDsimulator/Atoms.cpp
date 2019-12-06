#include "Atoms.h"
#include <iostream>
#include "dataType.h"

// The constructor initializes the position and velocity vectors and the
// distance matrix to the right size
Atoms::Atoms(int natoms, double m)
	: pos(natoms, vector<double>(3, 0)),
	vel(natoms, vector<double>(3, 0)),
	distances(natoms, vector<double>(natoms, 0)),
	reducedBondMatrix(natoms, vector<int>(natoms, 0)),
	bondTypes(0)
{
	// Initialize the number of atoms and the cell size
	nAtoms = natoms;
	cellLength = 0.0;
	mass = m;
	apm = natoms;
}

// The destructor deletes the memory of the position and velocity vectors
Atoms::~Atoms() {
	vector<vector<double>>().swap(pos);
	vector<vector<double>>().swap(vel);
	vector<vector<double>>().swap(distances);
	vector<vector<int>>().swap(reducedBondMatrix);
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

// Simple getter for apm
int Atoms::getApm() {
	return apm;
}

// Getter for the number of molecules in the simulation (or the number
// of repeated units)
int Atoms::getNM() {
	return nAtoms / apm;
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

// Only works if the atoms are not bonded in between the subcells
bool Atoms::isBonded(int i, int j) {
	int col_i = static_cast<int>(i / apm);
	int col_j = static_cast<int>(j / apm);
	if (col_i != col_j) {
		return false;
	}
	int red_i = i % apm;
	int red_j = j % apm;
	if (reducedBondMatrix[red_i][red_j] != 0) {
		return true;
	}
	return false;
}

double Atoms::getBondEnergy(int i, int j) {
	if (!isBonded(i, j)) return 0.0;
	bondT b = getBond(i, j);
	return b.getEnergy(distances[i][j]);
}

vector<double> Atoms::getBondForce(int i, int j) {
	if (!isBonded(i, j)) return { 0.0, 0.0, 0.0 };
	bondT b = getBond(i, j);
	return b.getForce(distances[i][j], pos[i], pos[j]);
}

bool Atoms::hasChangedPositions() {
	return positionsChanged;
}

// Casts the indexes into the reduced matrix
bondT Atoms::getBond(int i, int j) {
	int red_i = i % apm;
	int red_j = j % apm;
	int mIdx = reducedBondMatrix[red_i][red_j] + 1;
	return bondTypes[mIdx];
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

void Atoms::setBonds(vector<int> bonds, vector<double> ks, vector<double> r_es) {
	bondTypes.clear();  // Remove all existing bonds
	int counter = 1;  // set the first bond type to correspond to 1
	for (int i = 0; i < ks.size(); i++) {
		// Create a new bond and check whether its type already exists
		bondT b;
		b.ctor(ks[i], r_es[i]);
		bool alreadyAdded = false;
		int type;
		if (bondTypes.size() > 0) {
			for (int j = 0; j < bondTypes.size(); j++) {
				if (bondTypes[i] == b) {
					alreadyAdded = true;
					type = i + 1;
				}
			}
		}
		// If the type is new, then add it to the list, and increase the counter
		if (!alreadyAdded) {
			counter++;
			bondTypes.push_back(b);
			type = counter;
		}
		int second = i + 1;
		reducedBondMatrix[bonds[i]][bonds[second]] = type;
		reducedBondMatrix[bonds[second]][bonds[i]] = type;
	}
}

// Repeat the atoms object in the 3 cartesian direction, resulting in a
// NxNxN times bigger object.
void Atoms::repeat(int N) {
	if (N == 1) return;
	int size = nAtoms;
	int new_nMols = size / apm * N * N * N;
	resize(size / apm * N * N * N);
	int index = size;
	for (int x = 0; x < N; x++) {
		for (int y = 0; y < N; y++) {
			for (int z = 0; z < N; z++) {
				if (x + y + z == 0) continue;
				for (int i = 0; i < size; i++) {
					vector<double> p = pos[i];
					vector<double> t = {x * cellLength,
						y * cellLength, z * cellLength };
					for (int k = 0; k < 3; k++) {
						p[k] += t[k];
					}
					pos[index] = p;
					index++;
				}
			}
		}
	}
}

// Resizes the Atoms object and makes sure everything affected is updated
void Atoms::resize(int nMols) {
	int new_nAtoms = apm * nMols;
	nAtoms = new_nAtoms;
	pos.resize(new_nAtoms, vector<double>(3, 0));
	vel.resize(new_nAtoms, vector<double>(3, 0));
	distances.clear();
	distances.resize(new_nAtoms, vector<double>(new_nAtoms, 0));
	positionsChanged = true;
}