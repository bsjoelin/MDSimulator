#include "Ensemble.h"

Ensemble::Ensemble(Atoms* a, PotType potT, IntegratorType intT) {
	atoms = a;
	switch (potT)
	{
	case PotType::LJ:
		Pot = new LJ(atoms);
		break;
	default:
		Pot = new LJ(atoms);
		break;
	}
	switch (intT)
	{
	case IntegratorType::VERLET:
		Inte = new Verlet(atoms);
		break;
	default:
		Inte = new Verlet(atoms);
		break;
	}
}

Ensemble::~Ensemble() {
	delete &Pot, &Inte;
}

void Ensemble::calculate(double* U, vector<vector<double>>* F) {
	*U = Pot->getEnergy();
	*F = Pot->getForces();
}

void NVE::update(vector<vector<double>>* F) {
	Inte->update(F);
}