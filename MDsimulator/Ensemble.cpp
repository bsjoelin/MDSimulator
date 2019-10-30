#include "Ensemble.h"

Ensemble::Ensemble(Atoms* a, dataT* d)
	: forces(a->getSize(), vector<double>(3, 0))
{
	atoms = a;
	switch (d->PT)
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

	switch (d->IT)
	{
	case InteType::VERLET:
		InteEngine = new Verlet(atoms, &forces, d->dt_s);
		break;
	case InteType::VELVERLET:
		InteEngine = new VelVerlet(atoms, d->T_s, d->dt_s, d->tau_s_s);
		break;
	default:
		InteEngine = new Verlet(atoms, &forces, d->dt_s);
		break;
	}
}

Ensemble::~Ensemble() {
	delete &Pot, & InteEngine;
}

double Ensemble::calculate() {
	Pot->calculateDistances();
	forces = Pot->getForces();
	return Pot->getEnergy();
}

vector<vector<double>> Ensemble::getForces() {
	return Pot->getForces();
}

void Ensemble::printForces() {
	vector<double> av = { 0.0, 0.0, 0.0 };
	for (vector<double> f : forces) {
		for (int i = 0; i < 3; i++) {
			av[i] += f[i];
			cout << f[i] << ", ";
		}
		cout << endl;
	}
	cout << "Average: ";
	for (double d : av) {
		cout << d / atoms->getSize() << ", ";
	}
	cout << endl;
}

NVE::NVE(Atoms* a, dataT* d)
	: Ensemble(a, d) {}

void NVE::update() {
	InteEngine->update(atoms, &forces, this);
}