#include "Integrator.h"

Verlet::Verlet(Atoms* a, vector<vector<double>>* F, double diff_t)
	: oldPos(a->getSize(), vector<double>(3, 0))
{
	dt = diff_t;
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> q = vector<double>(3, 0);
		for (int j = 0; j < 3; j++) {
			double acc = (*F)[i][j];
			oldPos[i][j] = a->getPos(i)[j] - a->getVel(i)[j] * dt
				+ 1.0 / 2.0 * acc * dt * dt;
			q[j] = advancePos(a->getPos(i)[j], oldPos[i][j], acc);
		}
		a->setPos(i, q);
	}
}

Verlet::~Verlet() {
	oldPos.clear();
}

// The update leaves the positions one time step ahead of the velocities.
void Verlet::update(Atoms* a, vector<vector<double>>* F) {
	for (int i = 0; i < a->getSize(); i++) {
		vector<double> q = vector<double>(3, 0);
		vector<double> v = vector<double>(3, 0);
		for (int j = 0; j < 3; j++) {
			q[j] = advancePos(a->getPos(i)[j], oldPos[i][j], (*F)[i][j]);
			v[j] = advanceVel(q[j], oldPos[i][j]);
		}
		oldPos[i] = a->getPos(i);
		a->setPos(i, q);  // at t + dt
		a->setVel(i, v);  // at t
	}
}

double Verlet::advancePos(double q, double oldq, double F) {
	return 2.0 * q - oldq + F * dt * dt;
}

double Verlet::advanceVel(double newq, double oldq) {
	return 1.0 / (2.0 * dt) * (newq - oldq);
}


LeapFrog::LeapFrog(Atoms* a) {

}

void LeapFrog::update(vector<vector<double>>* F) {

}