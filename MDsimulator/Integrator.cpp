#include "Integrator.h"
#include "Ensemble.h"

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

Verlet::~Verlet() {
	oldPos.clear();
	nextPos.clear();
}

// The update leaves the positions and the velocities at the same time step.
void Verlet::update(Atoms* a, vector<vector<double>>* F, Ensemble* ens) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> nextnextq = vector<double>(3, 0);
		vector<double> v = vector<double>(3, 0);
		for (int j = 0; j < 3; j++) {
			nextnextq[j] = advancePos(nextPos[i][j], a->getPos(i)[j], (*F)[i][j]);
			v[j] = advanceVel(nextPos[i][j], oldPos[i][j]);
		}

		// Update all stored variables
		oldPos[i] = a->getPos(i);
		a->setPos(i, nextPos[i]);  // at t
		nextPos[i] = nextnextq;
		a->setVel(i, v);  // at t
	}
}

double Verlet::advancePos(double q, double oldq, double acc) {
	return 2.0 * q - oldq + acc * dt * dt;
}

double Verlet::advanceVel(double nextq, double oldq) {
	return (nextq - oldq) / (2.0 * dt);
}


VelVerlet::VelVerlet(Atoms* a, double temperature, double diff_t, double rel_t)
	: acc(a->getSize(), vector<double>(3, 0))
{
	dt = diff_t;
	T = temperature;
	Ms = 3.0 * a->getSize() * temperature * rel_t * rel_t;  // rel_t unitless
}

VelVerlet::~VelVerlet() {
	acc.clear();
}

void VelVerlet::update(Atoms* a, vector<vector<double>>* F, Ensemble* ens) {
	calculateAcceleration(a, F);
	if (Ms != 0.0) {  // if we are not using NVT, we just don't update zeta
		updateZeta(a, F);
	}
	updatePos(a);
	// get the forces from the advanced positions
	vector<vector<double>> nextForces = ens->getForces();
	updateVel(a, nextForces);
}

void VelVerlet::calculateAcceleration(Atoms* a, vector<vector<double>>* F) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++) {
			acc[i][j] = (*F)[i][j] - zeta * v[j];
		}
	}
}

void VelVerlet::updateZeta(Atoms* a, vector<vector<double>>* F) {
	// Calculate the sum of velocity times force
	double forcepos = 0;
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++)	{
			forcepos += v[i] * acc[i][j];
		}
	}

	// Update zeta
	zeta += dt / Ms * (a->getEnergy() - 3.0 * a->getSize() * T + dt * forcepos);
}

void VelVerlet::updatePos(Atoms* a) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> nextq = { 0.0, 0.0, 0.0 };
		vector<double> q = a->getPos(i);
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++)	{
			nextq[j] = q[j] + v[j] * dt + 1.0 / 2.0 * acc[i][j] * dt * dt;
		}
		a->setPos(i, nextq);
	}
}

void VelVerlet::updateVel(Atoms* a, vector<vector<double>> nF) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> nextv = { 0.0, 0.0, 0.0 };
		vector<double> v = a->getVel(i);
		for (int j = 0; j < 3; j++) {
			nextv[j] = v[j] + dt / 2.0 * (acc[i][j] + nF[i][j])
				/ (1.0 + zeta * dt / 2.0);
		}
		a->setVel(i, nextv);
	}
}

