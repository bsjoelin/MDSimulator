#ifndef _ensemble_h
#define _ensemble_h

#include "Potential.h"
#include "Integrator.h"
#include "dataType.h"

enum class EnsType { NVE, NVT };

// Class representing any ensemble (NVE, NVT, ...), which works as an interface
// with a few functions, which it passes on to its children.
class Ensemble
{
public:
	// Constructor and destrutor, which creates the correct Potential engine and
	// Integrator scheme, and takes care of destroying them when the ensemble is.
	Ensemble(Atoms* atoms, dataT* data);
	virtual ~Ensemble();

	// Public function for calculating energy and forces of the inherent Atoms object
	double calculate();
	// Calculate the pressure of the system
	double getPressure();
	// Public function for getting the forces from the Potential
	vector<vector<double>> getForces();
	// Print the forces vector to std::out
	void printForces();

	// Abstract function for updating the positions and velocities of the Atoms
	// object to the next time step. Must be implemented by children
	virtual double update() = 0;

protected:
	// The forces vector
	vector<vector<double>> forces;
	// Pointers to the inherent Atoms, Potential and Integrator objects
	Atoms* atoms;
	Potential* Pot;
	Integrator* InteEngine;
};

// Implementation of an NVE ensemble, which inherits from Ensemble
class NVE :
	public Ensemble
{
public:
	// Constructor
	NVE(Atoms* atoms, dataT* data);

	// Implementation of the abstrabt function update()
	double update();
};

// Implementation of an NVT ensemble, which inherits from Ensemble
class NVT :
	public Ensemble
{
public:
	// Constructor
	NVT(Atoms* atoms, dataT* data);

	// Implementation of the abstract function update()
	// Updates the positions and velocities and returns the energy of the
	// extended system
	double update();

private:
	double ln_s = 0;	// the natural logrithm of the scaling factor
	double zeta = 0;	// the friction coefficient
	double Ms;			// the thermal mass (reduced)
	double T;			// the reduced temperature
};

#endif // !_ensemble_h
