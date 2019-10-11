#ifndef _velocitymanager_h
#define _velocitymanager_h

#include "Atoms.h"
#include <random>
#include "time.h"

class VelocityManager {
public:
	VelocityManager();
	~VelocityManager();
	void initializeVelocities(Atoms* atoms, double temperature);

private:
	default_random_engine e;
	uniform_real_distribution<double> d;
	static constexpr double kB = 1.30864e-23;
	static constexpr double unitMass = 1.6726219e-27;

	double gaussianN(double mass, double temperature);
};

#endif // !_velocitymanager_h
