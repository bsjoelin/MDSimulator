#ifndef _integrator_h
#define _integrator_h

#include "Atoms.h"

// The Integrator class needs to know about the Ensemble class, but an include
// breaks things (circular inclusion), so instead we forward declare it
class Ensemble;

// Enumerator for all the implemented integrators
enum class InteType { VERLET, VELVERLET };

// Abstract class for an integrator for positions and velocities. Implementing
// classes must implement update()
class Integrator
{
public:
	// Abstract function for updating positions and velocities to next time step
	virtual void update(Atoms* atoms, vector<vector<double>>* forces,
		Ensemble* ens) = 0;

	// Functions for sending internal parameters up the chain for support for
	// a wide variety of ensembles
	void updateNvtParameters(double* ln_s, double* zeta);

protected:
	// The time step
	double dt = 0;

	// NVT parameters
	double T = 1;		// temperature
	double zeta = 0;	// 'friction' coefficient
	double Ms = 0.0;	// thermal mass (reduced)
	double ln_s = 0;	// natural log of scaling
};

// Implementation of the Verlet integrator scheme. Implements Integrator class
class Verlet :
	public Integrator
{
public:
	// Constructor and destructor
	Verlet(Atoms* atoms, vector<vector<double>>* forces, double dt);
	~Verlet();

	// Implementation of the abstract update() function
	void update(Atoms* atoms, vector<vector<double>>* forces, Ensemble* ens);

private:
	// The old and next positions have to be saved for the Verlet engine to have
	// the positions and velocities sync up
	vector<vector<double>> oldPos;
	vector<vector<double>> nextPos;

	// Functions for advancing a position and a velocity respectively
	double advancePos(double q, double oldq, double F);
	double advanceVel(double newq, double oldq);
};

// Implementation of the Velocity Verlet integrator scheme. This implements
// the Integrator class
class VelVerlet :
	public Integrator
{
public:
	// Constructor and destructor
	VelVerlet(Atoms* atoms, double T, double dt, double relaxation_time);
	virtual ~VelVerlet();

	// Implementation of the abstract update() function
	void update(Atoms* atoms, vector<vector<double>>* forces, Ensemble* ens);

private:
	vector<vector<double>> acc;  // saving the acceleration, since it's used often

	// Private functions for making the update work
	void calculateAcceleration(Atoms* atoms, vector<vector<double>>* forces);
	void updateZeta(Atoms* atoms, vector<vector<double>>* forces);
	void updatePos(Atoms* atoms);
	void updateVel(Atoms* atoms, vector<vector<double>> nextForces);
};

#endif // !_integrator_h
