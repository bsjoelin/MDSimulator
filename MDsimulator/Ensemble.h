#ifndef _ensemble_h
#define _ensemble_h

#include "Potential.h"
#include "Integrator.h"

class Ensemble
{
public:
	Ensemble(Atoms* atoms, PotType potT, InteType intT, double diff_t);
	~Ensemble();

	double calculate();
	virtual void update() = 0;


protected:
	Atoms* atoms;
	vector<vector<double>> forces;
	Potential* Pot;
	Integrator* Inte;
};

class NVE :
	public Ensemble
{
public:
	NVE(Atoms* a, PotType potT, InteType intT, double diff_t);

	void update();
};

#endif // !_ensemble_h
