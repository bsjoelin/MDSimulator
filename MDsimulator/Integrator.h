#ifndef _integrator_h
#define _integrator_h

#include "Atoms.h"

class Ensemble;

enum class InteType { VERLET, VELVERLET };

class Integrator
{
public:
	virtual void update(Atoms* atoms, vector<vector<double>>* forces, Ensemble* ens) = 0;

protected:
	double dt = 0;
};


class Verlet :
	public Integrator
{
public:
	Verlet(Atoms* atoms, vector<vector<double>>* forces, double dt);
	~Verlet();

	void update(Atoms* atoms, vector<vector<double>>* forces, Ensemble* ens);

private:
	vector<vector<double>> oldPos;
	vector<vector<double>> nextPos;

	double advancePos(double q, double oldq, double F);
	double advanceVel(double newq, double oldq);
};


class VelVerlet :
	public Integrator
{
public:
	VelVerlet(Atoms* atoms, double T, double dt, double relaxation_time);
	~VelVerlet();

	void update(Atoms* atoms, vector<vector<double>>* forces, Ensemble* ens);

private:
	double T;
	double zeta = 0;
	double Ms;
	vector<vector<double>> acc;

	void calculateAcceleration(Atoms* atoms, vector<vector<double>>* forces);
	void updateZeta(Atoms* atoms, vector<vector<double>>* forces);
	void updatePos(Atoms* atoms);
	void updateVel(Atoms* atoms, vector<vector<double>> nextForces);
};

#endif // !_integrator_h
