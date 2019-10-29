#include "Integrator.h"

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
}

// The update leaves the positions and the velocities at the same time step.
void Verlet::update(Atoms* a, vector<vector<double>>* F) {
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


LeapFrog::LeapFrog(Atoms* a) {

}

void LeapFrog::update(vector<vector<double>>* F) {

}