#ifndef _ensemble_h
#define _ensemble_h

#include "Potential.h"
#include "Integrator.h"
#include "dataType.h"


class Ensemble
{
public:
	Ensemble(Atoms* atoms, dataT* data);
	~Ensemble();

	double calculate();
	vector<vector<double>> getForces();
	void printForces();

	virtual void update() = 0;

protected:
	Atoms* atoms;
	vector<vector<double>> forces;
	Potential* Pot;
	Integrator* InteEngine;
};

class NVE :
	public Ensemble
{
public:
	NVE(Atoms* a, dataT* data);

	void update();
};

#endif // !_ensemble_h
