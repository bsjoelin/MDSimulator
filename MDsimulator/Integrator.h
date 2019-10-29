#ifndef _integrator_h
#define _integrator_h

#include "Atoms.h"

enum class InteType { VERLET };

class Integrator
{
public:
	virtual void update(Atoms* atoms, vector<vector<double>>* forces) = 0;
};


class Verlet :
	public Integrator
{
public:
	Verlet(Atoms* atoms, vector<vector<double>>* forces, double dt);
	~Verlet();

	void update(Atoms* atoms, vector<vector<double>>* forces);

private:
	vector<vector<double>> oldPos;
	vector<vector<double>> nextPos;
	double dt;

	double advancePos(double q, double oldq, double F);
	double advanceVel(double newq, double oldq);

};


class LeapFrog :
	public Integrator
{
public:
	LeapFrog(Atoms* atoms);

	void update(vector<vector<double>>* forces);
};

#endif // !_integrator_h
