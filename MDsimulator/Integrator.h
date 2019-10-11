#ifndef _integrator_h
#define _integrator_h

#include "Atoms.h"

enum class IntegratorType { VERLET };

class Integrator
{
public:
	virtual void update(vector<vector<double>>* forces) = 0;
};

class Verlet :
	public Integrator
{
public:
	Verlet(Atoms* atoms);

	void update(vector<vector<double>>* forces);
};

#endif // !_integrator_h
