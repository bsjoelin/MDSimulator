#include "Ensemble.h"
#include <iostream>

// Constructor for any Ensemble, which assigns the Atoms object and creates
// the wanted Potential and Integrator objects with the needed parameters.
Ensemble::Ensemble(Atoms* a, dataT* d)
	: forces(a->getSize(), vector<double>(3, 0))  // initialize forces vector
{
	atoms = a;  // Assign Atoms pointer
	// Switch on the Potential type, and create the proper one
	switch (d->PT)
	{
	case PotType::LJ:
		Pot = new LJ(atoms, d->rhoN, d->r_co);
		break;
	// default is a Lennard-Jones Potential
	default:
		Pot = new LJ(atoms, d->rhoN, d->r_co);
		break;
	}

	// Get the forces from the Potential
	forces = Pot->getForces();

	// Switch on the Integrator type and create the proper one
	switch (d->IT)
	{
	case InteType::VERLET:
		InteEngine = new Verlet(atoms, &forces, d->dt_s);
		break;
	case InteType::VELVERLET:
		InteEngine = new VelVerlet(atoms, d->T_s, d->dt_s, d->tau_s_s);
		break;
	// Default is the Verlet, which is only really for NVE
	default:
		InteEngine = new Verlet(atoms, &forces, d->dt_s);
		break;
	}
}

// Destructor that deletes the forces vector and the Potential and Integrator
// objects
Ensemble::~Ensemble() {
	vector<vector<double>>().swap(forces);
	delete &Pot, &InteEngine;
}

Ensemble* Ensemble::createEnsemble(Atoms* a, dataT* d) {
	switch (d->ET)
	{
	case EnsType::NVE:
		return new NVE(a, d);
	case EnsType::NVT:
		return new NVT(a, d);
	default:
		return new NVE(a, d);
	}
}

// Ask the Potential to calculate the distances between atoms, so the potential
// energy and the forces can be calculated
double Ensemble::calculate() {
	forces = Pot->getForces();
	return Pot->getEnergy();
}

double Ensemble::getPressure() {
	return (2 * atoms->getEnergy() + Pot->getSumForcesInteraction())
		/ (3 * pow(atoms->getCellLength(), 3.0))
		+ Pot->getPressureCorrection();
}

// Wrapper for getting the forces from the potential, when the stored forces
// are not the ones needed
vector<vector<double>> Ensemble::getForces() {
	return Pot->getForces();
}

// Simple printing function for printing the forces to the console
void Ensemble::printForces() {
	vector<double> av = { 0.0, 0.0, 0.0 };
	for (vector<double> f : forces) {
		for (int i = 0; i < 3; i++) {
			av[i] += f[i];
			cout << f[i] << ", ";
		}
		cout << endl;
	}
	// Prints average of forces as well
	cout << "Average: ";
	for (double d : av) {
		cout << d / atoms->getSize() << ", ";
	}
	cout << endl;
}

// The constructor for NVE calls the Ensemble constructor
NVE::NVE(Atoms* a, dataT* d)
	: Ensemble(a, d) {}

// The update() function asks the Integrator to update
double NVE::update() {
	InteEngine->update(atoms, &forces, this);
	return 0;  // There is no extended system for NVE, so return zero
}

// The constructor for NVT calls the Ensemble constructor
NVT::NVT(Atoms* a, dataT* d)
	: Ensemble(a, d) 
{
	T = d->T_s;
	Ms = 3.0 * a->getSize() * T * d->tau_s_s * d->tau_s_s;  // rel_t unitless
}

// The update() function asks the Integrator to update, and the returns the
// energy of the extended system
double NVT::update() {
	InteEngine->update(atoms, &forces, this);
	// add the energy from the extended system
	InteEngine->updateNvtParameters(&ln_s, &zeta);
	return zeta * zeta * Ms / 2.0 + 3.0 * atoms->getSize() * T * ln_s;
}
