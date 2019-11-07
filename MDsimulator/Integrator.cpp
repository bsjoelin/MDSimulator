#include "Integrator.h"
#include "Ensemble.h"


void Integrator::updateNvtParameters(double* _ln_s, double* _zeta) {
	*_ln_s = ln_s;
	*_zeta = zeta;
}


// The constructor initializes and populates the new and old positions vectors
Verlet::Verlet(Atoms* a, vector<vector<double>>* F, double diff_t)
	: oldPos(a->getSize(), vector<double>(3, 0)),
	nextPos(a->getSize(), vector<double>(3, 0))
{
	dt = diff_t;
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> q = vector<double>(3, 0);
		for (int j = 0; j < 3; j++) {
			double acc = (*F)[i][j];
			oldPos[i][j] = a->getPos(i)[j] - a->getVel(i)[j] * dt
				+ 1.0 / 2.0 * acc * dt * dt;
			nextPos[i][j] = advancePos(a->getPos(i)[j], oldPos[i][j], acc);
		}
	}
}

// The destructor releases memory from the internal vectors
Verlet::~Verlet() {
	vector<vector<double>>().swap(oldPos);
	vector<vector<double>>().swap(nextPos);
}

// The update leaves the positions and the velocities at the same time step.
void Verlet::update(Atoms* a, vector<vector<double>>* F, Ensemble* ens) {
	// Run through all the atoms and update their positions and velocities
	for (int i = 0; i < a->getSize(); i++) {
		// Positions at t + 2dt
		vector<double> nextnextq = vector<double>(3, 0);
		// Empty velocity vector
		vector<double> v = vector<double>(3, 0);
		for (int j = 0; j < 3; j++) {
			nextnextq[j] = advancePos(nextPos[i][j], a->getPos(i)[j], (*F)[i][j]);
			v[j] = advanceVel(nextPos[i][j], oldPos[i][j]);
		}

		// Update all stored variables
		oldPos[i] = a->getPos(i);  // q(t - dt) = q(t)
		a->setPos(i, nextPos[i]);  // q(t) = q(t + dt)
		nextPos[i] = nextnextq;    // q(t + dt) = q(t + 2dt)
		a->setVel(i, v);
	}
}

double Verlet::advancePos(double q, double oldq, double acc) {
	// q(t + dt) = 2q(t) - q(t - dt) + a(t) * dt * dt
	return 2.0 * q - oldq + acc * dt * dt;
}

double Verlet::advanceVel(double nextq, double oldq) {
	// v(t + dt) = (q(t + dt) - q(t - dt)) / 2dt
	return (nextq - oldq) / (2.0 * dt);
}


// The constructor calculates the thermal mass and initializes internal
// memory members
VelVerlet::VelVerlet(Atoms* a, double temperature, double diff_t, double rel_t)
	: acc(a->getSize(), vector<double>(3, 0))
{
	dt = diff_t;
	T = temperature;
	Ms = 3.0 * a->getSize() * temperature * rel_t * rel_t;  // rel_t unitless
}

// The destructor releases the memory of the acceleration vector
VelVerlet::~VelVerlet() {
	vector<vector<double>>().swap(acc);
}

// The update() function works for both NVE and NVT, i.e. in NVT reduecs to
// NVE when zeta = 0. So regardsless of the value of zeta, the function updates
// the positions and velocities to the next time step
void VelVerlet::update(Atoms* a, vector<vector<double>>* F, Ensemble* ens) {
	// Calculate all the accelarations
	calculateAcceleration(a, F);
	if (Ms != 0.0) {  // if we are not using NVT, we just don't update zeta
		updateZeta(a, F);
	}
	// Update the postions in the Atoms object
	updatePos(a);
	// The Velocity Verlet method use the forces from the next iteration, so
	// we recalculate the forces from the now updated positions
	vector<vector<double>> nextForces = ens->getForces();
	updateVel(a, nextForces);
}

// Function for calculating all the accelerations
void VelVerlet::calculateAcceleration(Atoms* a, vector<vector<double>>* F) {
	// In the case of NVE, zeta = 0, so the calculation reduces to a = F / m
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++) {
			acc[i][j] = (*F)[i][j] - zeta * v[j];
		}
	}
}

// Function for updating the friction coefficient
void VelVerlet::updateZeta(Atoms* a, vector<vector<double>>* F) {
	// Calculate the sum of velocity times acceleration
	double forcepos = 0;
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++)	{
			forcepos += v[j] * acc[i][j];
		}
	}
	// Calculate the rate of change of zeta
	double dotZeta = (2.0 * a->getEnergy() - 3.0 * a->getSize() * T) / Ms;

	// Update ln(s)
	ln_s += zeta * dt + 1.0 / 2.0 * dotZeta * dt * dt;

	// Update zeta
	zeta += dt * (dotZeta + dt / Ms * forcepos);
}

// Function for updating the positions
void VelVerlet::updatePos(Atoms* a) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> nextq = { 0.0, 0.0, 0.0 };
		vector<double> q = a->getPos(i);
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++)	{
			// q(t + dt) = q(t) + v(t) * dt + 1 / 2 * a(t) * dt * dt
			nextq[j] = q[j] + v[j] * dt + 1.0 / 2.0 * acc[i][j] * dt * dt;
		}
		a->setPos(i, nextq);
	}
}

// Function for updating the velocities
void VelVerlet::updateVel(Atoms* a, vector<vector<double>> nF) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> nextv = { 0.0, 0.0, 0.0 };
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++) {
			// v(t + dt) = (v(t) + 0.5 * dt * (a(t) + a(t + dt)) 
			//    / (1 + zeta(t + dt) * 0.5 * dt
			nextv[j] = (v[j] + dt / 2.0 * (acc[i][j] + nF[i][j]))
				/ (1.0 + zeta * dt / 2.0);
		}
		a->setVel(i, nextv);
	}
}
