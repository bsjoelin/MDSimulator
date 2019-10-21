#ifndef _potential_h
#define _potential_h

#include "Atoms.h"

enum class PotType { LJ };

class Potential
{
public:
	
	Potential(Atoms* atoms);

	virtual double getEnergy() = 0;
	virtual vector<vector<double>> getForces() = 0;
	void calculateDistances();

protected:
	Atoms* atoms;
	vector<vector<double>> dist;
};


class LJ :
	public Potential
{
public:
	LJ(Atoms* a);
	~LJ();

	double getEnergy();
	vector<vector<double>> getForces();
	void printDistances();
	void printForces(vector<vector<double>> F);
};

#endif // !_potential_h
