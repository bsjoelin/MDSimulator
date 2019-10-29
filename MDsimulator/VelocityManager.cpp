#define _USE_MATH_DEFINES

#include "VelocityManager.h"

VelocityManager::VelocityManager() {}

VelocityManager::~VelocityManager() {}

void VelocityManager::initializeVelocities(Atoms* atoms, double T) {
	vector<double> v = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < atoms->getSize(); i++) {
		for (int j = 0; j < 3; j++) {
			v[j] = gaussianN(atoms->mass, T);
		}
		atoms->setVel(i, v);
	}
	atoms->centerVel();

	double T_calc = 2 / (3.0 * atoms->getSize()) * atoms->getEnergy();
	double correction = pow(T / T_calc, 0.5);

	for (int i = 0; i < atoms->getSize(); i++) {
		vector<double> v_c = atoms->getVel(i);
		for (int j = 0; j < 3; j++) {
			v_c[j] *= correction;
		}
		atoms->setVel(i, v_c);
	}
}

double VelocityManager::gaussianN(double m, double T) {
	// The gaussian variance is calculated
	double variance = T;  // * m[i];
	// The two polar coordinates for a gaussian distribution
	double phi = 2 * M_PI * d(e);
	double s = -variance * log(1.0 - d(e));
	// Return the cartesian gaussian
	return pow(2 * s, 0.5) * cos(phi);
}
