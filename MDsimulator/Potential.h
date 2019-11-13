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
	Potential(Atoms* atoms, double numberDensity, double radialCutOff);
	virtual ~Potential();

	// Function for getting the potential energy of the collection of atoms
	virtual double getEnergy() = 0;
	// Function for getting the forces on every atom
	virtual vector<vector<double>> getForces() = 0;
	// Calculate the pressure tail correction resulting from the cut-off
	virtual double getPressureCorrection() = 0;

	// Function for retrieving the sum of force interactions ((r_i - r_j) * F_ji)
	double getSumForcesInteraction();

// The following menbers are protected, so they are inherited by implementing
// classes
protected:
	Atoms* atoms;
	double numberDensity;	// number density (dimension-less)
	double r_c;				// radial cut-off (dimension-less)
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
	LJ(Atoms* a, double numberDensity, double radialCutoff);

	// Implements the abstract functions getEnergy() and getForces()
	double getEnergy();
	vector<vector<double>> getForces();
	// Calculate the pressure tail correction resulting from the cut-off
	double getPressureCorrection();
	
	// Helper functions for writing the force vector to console
	void printForces(vector<vector<double>> F);

private:
	double cutoffEnergy;	// the energy at the cut-off
	double cutoffForce;		// the force at the cut-off

	// Calculate the energy between a single pair, and handle cut-off
	double calculateEnergy(double distance);
	// Calculate the energy tail correction resulting from the cut-off
	double calculateEnergyCorrection();
};

#endif // !_potential_h
