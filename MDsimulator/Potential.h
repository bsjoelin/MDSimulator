#ifndef _potential_h
#define _potential_h

#include "Atoms.h"

// Enumerator containing the implemented potential types
enum class PotType { LJ };

// Abstract class for making a potential for atom interaction. Implementing
// classes must implement getEnergy() and getForces()
class Potential
{
public:
	// Constructor and destructor
	Potential(Atoms* atoms);
	virtual ~Potential();

	// Function for getting the potential energy of the collection of atoms
	virtual double getEnergy() = 0;
	// Function for getting the forces on every atom
	virtual vector<vector<double>> getForces() = 0;

	// Function for retrieving the sum of force interactions ((r_i - r_j) * F_ji)
	double getSumForcesInteraction();
	// Let the potential recalculate distances and forces
	void reset();

// The following menbers are protected, so they are inherited by implementing
// classes
protected:
	Atoms* atoms;
	// This is for the allowing the recalculation of the forces
	bool forcesCalculated = false;
	// Keep the results of distances and forces in memory to reduce
	// computational cost
	double sumForceInteractions = 0;
	vector<vector<double>> dist;
	vector<vector<double>> forces;
};

// Implementation of the Potential class with a Lennard-Jones 12-6 potential
class LJ :
	public Potential
{
public:
	LJ(Atoms* a);

	// Implements the abstract functions getEnergy() and getForces()
	double getEnergy();
	vector<vector<double>> getForces();
	// Helper functions for writing the distance matrix or force vector to console
	void printDistances();
	void printForces(vector<vector<double>> F);
};

#endif // !_potential_h
