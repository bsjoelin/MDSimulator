#include "Potential.h"
#include <iostream>

Potential::Potential(Atoms* a) {
	atoms = a;
}

void Potential::getDistances(Atoms* atoms, vector<vector<double>>* dist) {
	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			vector<double> pos1 = atoms->getPos(i);
			vector<double> pos2 = atoms->getPos(j);
			double r = 0.0;
			for (int k = 0; k < 3; k++) {
				double diff = pos1[k] - pos2[k];
				r += diff * diff;
			}
			r = pow(r, 0.5);
			(*dist)[i][j] = r;
			(*dist)[j][i] = r;
		}
	}
}

LJ::LJ(Atoms* a) 
	: dist(a->getSize(), vector<double>(a->getSize(), 0)),
	Potential(a) 
{
	atoms = a;
}

LJ::~LJ() {
	delete& dist;
}


double LJ::getEnergy() {
	getDistances(atoms, &dist);
	printDistances();

	double U = 0.0;

	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			U += 4.0 * (pow(1.0 / dist[i][j], 12.0) - pow(1.0 / dist[i][j], 6.0));
		}
	}

	return U;
}

vector<vector<double>> LJ::getForces() {
	
	vector<vector<double>> F(atoms->getSize(), vector<double>(3, 0));
	for (int i = 0; i < atoms->getSize() - 1; i++) {
		for (int j = i + 1; j < atoms->getSize(); j++) {
			// scalar force divided by distance
			double pf = 48 / dist[i][j] * (pow(1 / dist[i][j], 14)
				- 0.5 * pow(1 / dist[i][j], 12));
			for (int k = 0; k < 3; k++)	{
				double F_jia = pf * (atoms->getPos(i)[k] - atoms->getPos(j)[k]);
				
				F[i][k] += F_jia;
				F[j][k] -= F_jia;
			}
		}
	}
	printForces(F);
	return F;
}

void LJ::printDistances() {
	for (vector<double> p : dist) {
		for (double d : p) {
			cout << d << ", ";
		}
		cout << "\n";
	}
}

void LJ::printForces(vector<vector<double>> F) {
	for (vector<double> f : F) {
		for (double a : f) {
			cout << a << ", ";
		}
		cout << endl;
	}
}