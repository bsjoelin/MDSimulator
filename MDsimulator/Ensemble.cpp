#include "Ensemble.h"

Ensemble::Ensemble(Atoms* a, PotType potT, InteType intT, double diff_t)
	: forces(a->getSize(), vector<double>(3, 0))
{
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
	forces = Pot->getForces();

	switch (intT)
	{
	case InteType::VERLET:
		Inte = new Verlet(atoms, &forces, diff_t);
		break;
	default:
		Inte = new Verlet(atoms, &forces, diff_t);
		break;
	}
}

Ensemble::~Ensemble() {
	delete &Pot, &Inte;
}

double Ensemble::calculate() {
	Pot->calculateDistances();
	forces = Pot->getForces();
	return Pot->getEnergy();
}

NVE::NVE(Atoms* a, PotType potT, InteType intT, double diff_t)
	: Ensemble(a, potT, intT, diff_t) {}

void NVE::update() {
	Inte->update(atoms, &forces);
}