#include "Atoms.h"


Atoms::Atoms(int natoms)
	: pos(natoms, vector<double>(3, 0)),
	vel(natoms, vector<double>(3, 0))
{
	nAtoms = natoms;
	cellLength = 0.0;
}

Atoms::~Atoms() {
	pos.clear();
	vel.clear();
}

void Atoms::print(bool printAverage) {
	vector<double> R = { 0.0, 0.0, 0.0 };
	for (vector<double> p : pos) {
		for (int i = 0; i < 3; i++) {
			R[i] += p[i];
			cout << p[i] << ", ";
		}
		cout << "\n";
	}
	if (printAverage) {
		cout << "Average: ";
		for (double d : R) {
			cout << d / nAtoms << ", ";
		}
		cout << endl;
	}
}

void Atoms::printVel() {
	vector<double> av = { 0.0, 0.0, 0.0 };
	for (vector<double> v : vel) {
		for (int i = 0; i < 3; i++) {
			av[i] += v[i];
			cout << v[i] << ", ";
		}
		cout << "\n";
	}
	cout << "Average: ";
	for (double d : av) {
		cout << d / nAtoms << ", ";
	}
	cout << endl;
}

void Atoms::center() {
	vector<double> R = { 0.0, 0.0, 0.0 };
	double M = 0;
	// find the center of mass (COM)
	for (vector<double> p : pos) {
		for (int i = 0; i < 3; i++) {
			R[i] += mass * p[i];
		}
		M += mass;
	}
	for (int i = 0; i < 3; i++)
	{
		R[i] /= M;
	}

	// adjust all positions, so they lie around COM
	for (vector<double> &p : pos) {
		for (int i = 0; i < 3; i++)	{
			p[i] -= R[i];
		}
	}
}

void Atoms::centerVel() {
	vector<double> av = { 0.0, 0.0, 0.0 };
	double M = 0.0;
	for (vector<double> v : vel) {
		for (int i = 0; i < 3; i++)	{
			av[i] += mass * v[i];
		}
		M += mass;
	}
	for (int i = 0; i < 3; i++) {
		av[i] /= M;
	}

	for (vector<double> &v : vel) {
		for (int i = 0; i < 3; i++)	{
			v[i] -= av[i];
		}
	}
}

int Atoms::getSize() {
	return nAtoms;
}

double Atoms::getCellLength() {
	return cellLength;
}

vector<double> Atoms::getPos(int i) {
	return pos[i];
}

vector<double> Atoms::getVel(int i) {
	return vel[i];
}

double Atoms::getEnergy() {
	double K = 0.0;
	for (vector<double> v : vel) {
		for (double vj : v) {
			K += vj * vj;
		}
	}
	return 1.0 / 2.0 * K;
}

void Atoms::setPos(int i, vector<double> r) {
	pos[i] = r;
}

void Atoms::setVel(int i, vector<double> r) {
	vel[i] = r;
}

void Atoms::setCellLength(double length) {
	if (length > 0.0) {
		cellLength = length;
	}
}