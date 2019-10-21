#define _USE_MATH_DEFINES

#include "VelocityManager.h"
#include <cmath>
#include <iostream>

VelocityManager::VelocityManager(){
	
}

VelocityManager::~VelocityManager() {

}

void VelocityManager::initializeVelocities(Atoms* atoms, double T) {
	vector<double> v = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < atoms->getSize(); i++) {
		for (int j = 0; j < 3; j++) {
			v[j] = gaussianN(atoms->mass, T);
		}
		atoms->setVel(i, v);
	}
	atoms->centerVel();
	//atoms->printVel();
}

double VelocityManager::gaussianN(double m, double T) {
	// The gaussian variance is calculated
	double variance = kB * T / (m * unitMass);
	// The two polar coordinates for a gaussian distribution
	double phi = 2 * M_PI * d(e);
	double s = -variance * log(1.0 - d(e));
	// Return the cartesian gaussian
	return pow(2 * s, 0.5) * cos(phi);
}
