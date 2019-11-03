#ifndef _velocitymanager_h
#define _velocitymanager_h

#include "Atoms.h"

// Static class for initializing the velocities to a proper gaussian distibution
// and make sure that it corresponds to the given temperature.
class VelocityManager {
public:
	// Static function for initializing the velocities of an Atoms object to
	// correspond to a proper gaussian distribution, which has a kinetic energy
	// in accordance to the temperature given.
	static void initializeVelocities(Atoms* atoms, double temperature);

private:
	// Private helper function for generating a gaussian number from a uniform
	// distribution in [0, 1].
	static double gaussianN(double mass, double temperature);
};

#endif // !_velocitymanager_h
