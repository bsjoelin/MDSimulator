#include "Ensemble.h"

Ensemble::Ensemble(Atoms* a, PotType potT, InteType intT, double diff_t) {
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

	Pot->calculateDistances();
	vector<vector<double>> F = Pot->getForces();

	switch (intT)
	{
	case InteType::VERLET:
		Inte = new Verlet(atoms, &F, diff_t);
		break;
	default:
		Inte = new Verlet(atoms, &F, diff_t);
		break;
	}
}

Ensemble::~Ensemble() {
	delete &Pot, &Inte;
}

void Ensemble::calculate(double* U, vector<vector<double>>* F) {
	Pot->calculateDistances();
	*U = Pot->getEnergy();
	*F = Pot->getForces();
}

NVE::NVE(Atoms* a, PotType potT, InteType intT, double diff_t)
	: Ensemble(a, potT, intT, diff_t) {}

void NVE::update(vector<vector<double>>* F) {
	Inte->update(atoms, F);
}