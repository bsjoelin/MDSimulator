#include "Atoms.h"
#include <iostream>


Atoms::Atoms(int natoms) {
	nAtoms = natoms;
	cellLength = 0.0;
	pos.resize(natoms, vector<double>(3, 0));
	vel.resize(natoms, vector<double>(3, 0));
}

Atoms::~Atoms() {
	pos.clear();
	oldPos.clear();
	vel.clear();
	oldVel.clear();
}

void Atoms::print() {
	vector<double> R = { 0.0, 0.0, 0.0 };
	for (vector<double> p : pos) {
		for (int i = 0; i < 3; i++) {
			R[i] += p[i];
			cout << p[i] << ", ";
		}
		cout << "\n";
	}
	cout << "Average: ";
	for (double d : R) {
		cout << d << ", ";
	}
	cout << endl;
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
		cout << d << ", ";
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

void Atoms::resize(int n) {
	nAtoms = n;
	pos.resize(n, vector<double>(3, 0));
	vel.resize(n, vector<double>(3, 0));
}

int Atoms::getSize() {
	return nAtoms;
}

vector<double> Atoms::getPos(int i) {
	return pos[i];
}

vector<double> Atoms::getVel(int i) {
	return vel[i];
}

void Atoms::setPos(int i, vector<double> r) {
	pos[i] = r;
}

void Atoms::setVel(int i, vector<double> r) {
	vel[i] = r;
}

void Atoms::setCellLength(double length) {
	if (cellLength > 0.0) {
		cellLength = length;
	}
}
