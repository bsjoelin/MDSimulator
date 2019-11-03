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
		Pot = new LJ(atoms);
		break;
	// default is a Lennard-Jones Potential
	default:
		Pot = new LJ(atoms);
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

// Ask the Potential to calculate the distances between atoms, so the potential
// energy and the forces can be calculated
double Ensemble::calculate() {
	forces = Pot->getForces();
	return Pot->getEnergy();
}

double Ensemble::getPressure() {
	return (2 * atoms->getEnergy() + Pot->getSumForcesInteraction())
		/ (3 * pow(atoms->getCellLength(), 3.0));
}

// Wrapper for getting the forces from the potential, when the stored forces
// are not the ones needed
vector<vector<double>> Ensemble::getForces() {
	return Pot->getForces();
}

// Lets the Potential recalculate distances and forces
void Ensemble::resetPot() {
	Pot->reset();
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
	: Ensemble(a, d) {}

// The update() function asks the Integrator to update, and the returns the
// energy of the extended system
double NVT::update() {
	InteEngine->update(atoms, &forces, this);
	// add the energy from the extended system
	return 0;  // Must be the energy of the extended system
}