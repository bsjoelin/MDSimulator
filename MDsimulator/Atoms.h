#ifndef _atoms_h
#define _atoms_h

#include <vector>

using namespace std;

// A class that works as a container for single atoms, so can be a representation
// of any number of atoms and/or molecules. It keeps track of position, velocity
// and maybe mass at some point.
class Atoms {
public:
	// Constructor and destructor for object. The constructur takes the number of
	// atoms it will contain.
	Atoms(int natoms, double mass);
	virtual ~Atoms();

	// Functions for centering the positions and velocities, so the weighted
	// average is zero in all cartesian coordinates.
	void center();
	void centerVel();

	// Functions for printing the positions, velocities and distances to the console
	void print();
	void printVel();
	void printDistances();

	// Function for calculating all the interatomic distances
	vector<vector<double>> getDistances();

	// Getter functions for the object members
	int getSize();  // Get number of atoms
	double getCellLength();  // Get the side length of the cell [dimensionless]
	double getMass();  // Get the mass of the atoms
	vector<double> getPos(int i);  // Get the position vector of atom i
	vector<double> getVel(int i);  // Get the velocity vector of atom i
	double getEnergy();  // Get the kinetic energy of all the atoms
	bool hasChangedPositions();  // Get whether or not the positions have changed

	// Setter functions for the object members
	void setPos(int i, vector<double> r);  // Set the position vector of atom i as r
	void setVel(int i, vector<double> r);  // Set the velocity vector of atom i as r
	void setCellLength(double length);  // Set the side length of the cell

private:
	// Used for determining when to recalculate distance matrix
	bool positionsChanged = true;
	int nAtoms;  // The number of atoms
	double mass;  // The mass of the atoms
	double cellLength;  // The side length of the cell
	vector<vector<double> > pos, vel;  // Position and velocity vectors
	vector<vector<double>> distances;
};

#endif // !_atoms_h