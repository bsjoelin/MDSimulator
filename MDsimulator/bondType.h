#ifndef _bondtype_h
#define _bondtype_h
#include "math.h"
#include <vector>

struct bondT {
	double k_s = 0;		// Reduced bond force constant
	double r_eq_s = 1;	// Reduced bond equilibrium distance

	// Fake constructor
	void ctor(double k, double r) {
		k_s = k;
		r_eq_s = r;
	}

	// Get the bond energy
	double getEnergy(double dist) {
		return 0.5 * k_s * pow(dist - r_eq_s, 2.0);
	}

	// Get the 3D bond force vector.
	std::vector<double> getForce(double dist, std::vector<double> p1, std::vector<double> p2) {
		double force = k_s * (dist - r_eq_s) / dist;
		std::vector<double> v = { force, force, force };
		for (int k = 0; k < 3; k++) {
			double diff = p2[k] - p1[k];
			v[k] *= diff;
		}
		return v;
	}

	// Override for the 'equals' operator
	bool operator ==(bondT other) {
		if (k_s == other.k_s && r_eq_s == other.r_eq_s) {
			return true;
		}
		return false;
	}
};

#endif // !_bondtype_h

