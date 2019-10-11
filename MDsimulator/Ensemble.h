#ifndef _ensemble_h
#define _ensemble_h

#include "Potential.h"
#include "Integrator.h"

class Ensemble
{
public:
	Ensemble(Atoms* atoms, PotType potT, IntegratorType intT);
	~Ensemble();

	void calculate(double* energy, vector<vector<double>>* forces);
	virtual void update(vector<vector<double>>* forces) = 0;


protected:
	Atoms* atoms;
	Potential* Pot;
	Integrator* Inte;
};

class NVE :
	public Ensemble
{
public:
	void update(vector<vector<double>>* forces);
private:

};

#endif // !_ensemble_h
