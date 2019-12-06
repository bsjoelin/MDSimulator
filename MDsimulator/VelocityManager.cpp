#define _USE_MATH_DEFINES
#include "VelocityManager.h"
#include <random>
#include "time.h"

// A struct needed to be able to use the class in a static way, which can
// produce a random real number in the range [0, 1].
struct uniform_double {
	double random = generate();
	static double generate() {
		static default_random_engine e{
			static_cast<long unsigned int>(time(0))  // Pseudo random seed
		};
		static uniform_real_distribution<double> d{0, 1};
		return d(e);
	}
};

// The static initializeVelocities() function generates the velocities from
// a gaussian distribution and makes sure that the center of velocity is zero.
void VelocityManager::initializeVelocities(Atoms* atoms, double T) {
	vector<double> v = { 0.0, 0.0, 0.0 };
	// Loop through all velocities, and generate them
	for (int i = 0; i < atoms->getSize(); i++) {
		for (int j = 0; j < 3; j++) {
			v[j] = gaussianN(atoms->getMass(), T);
		}
		atoms->setVel(i, v);
	}
	// Center the velocity
	atoms->centerVel();

	// Calculate the actual (instantaneous) temperature
	double T_calc = 2 / (3.0 * atoms->getSize()) * atoms->getEnergy();
	// Calculate the correction/scaling factor to achieve desired temperature
	double correction = pow(T / T_calc, 0.5);

	// Scale all the velocites with this correction factor
	for (int i = 0; i < atoms->getSize(); i++) {
		vector<double> v_c = atoms->getVel(i);
		for (int j = 0; j < 3; j++) {
			v_c[j] *= correction;
		}
		atoms->setVel(i, v_c);
	}
}

// The gaussianN() function calculates a number from a gaussian distribution
// from a uniformly distributed random number in [0, 1].
double VelocityManager::gaussianN(double m, double T) {
	// The gaussian variance is calculated
	double variance = T;  // * m[i];
	// The two polar coordinates for a gaussian distribution
	uniform_double ra = uniform_double();
	double phi = 2 * M_PI * ra.random;
	double s = -variance * log(1.0 - ra.random);
	// Return the cartesian gaussian
	return pow(2 * s, 0.5) * cos(phi);
}
