#include "CellBuilder.h"
#include <iostream>
#include <cmath>
#include "math.h"

// The static buildCell() function determines the correct lattice sytem to
// build and sends the necessary parameters on to the proper functions
void CellBuilder::buildCell(Atoms* atoms, int nMols, double density) {
	int N = nMols;  // the number of molecules
	double n = pow(N, (1.0 / 3.0));  // the cube root of N
	// Calculate the cell length as described in the Atoms object
	double cellLength = pow(N / density, (1.0 / 3.0));

	// The following determines the correct system to build by calculating
	// the prefactor in the equation
	//		N = c * n^3
	// where c = 1 for a sc system, c = 2 for bcc, and c = 4 for fcc, and then
	// determines when (n/c)^(1/3) is an integer.
	// The use of modulos (modf) is to determine when (n/c)^(1/3) is an integer,
	// since modf(double d, double* i) returns the fractional part of d, and
	// places the integral part in i. So if modf(n * c^(1/3)) << 1, we have
	// found the correct system. The small added number is to ensure we don't
	// get modf(n * c^(1/3)) = 0.9999...

	double temp = 0.0;  // the integer
	if (abs(modf(n + 0.0001, &temp)) < 0.001) {
		BuildSCUnitCell(atoms, temp, cellLength);
	}
	else if (abs(modf(n * bccFactor + 0.001, &temp)) < 0.01) {
		buildBCCUnitCell(atoms, temp, cellLength);
	}
	else if (abs(modf(n * fccFactor + 0.001, &temp)) < 0.01) {
		buildFCCUnitCell(atoms, temp, cellLength);
	}
	// If c is neither 1, 2 or 4, we can assume, that the number is not a
	// magic number, and we just round N down to the nearest sc system.
	else {
		N = int(pow(floor(n), 3.0));
		cout << "Number of atoms is not a magic number! using "
			<< N << " instead." << endl;

		// Resize the atoms object and calculate the new cell length.
		double new_n = pow(N, (1.0 / 3.0));
		cellLength = pow(N / density, (1.0 / 3.0));  // new cell size
		BuildSCUnitCell(atoms, new_n, cellLength);
	}
	// Give the cell size to the Atoms object
	atoms->setCellLength(cellLength);
	// Center the atoms
	atoms->center();
}


// All the following functions work by:
// 1) Determine the number of subcells (n^3)
// 2) Calculate the half side length of a subcell
// 3) Place the atoms in the unit cell as the lattice type prescribes
// 4) set the cell length as the unit cell length

void CellBuilder::BuildSCUnitCell(Atoms* atoms, double N, double length) {
	int n = static_cast<int>(N);
	double dHalfCell = 0.5 * length / n;
	atoms->setCellLength(2.0 * dHalfCell);
	atoms->repeat(n);
	cout << "Chose SC!" << endl;
}

void CellBuilder::buildBCCUnitCell(Atoms* atoms, double N, double length) {
	int n = static_cast<int>(N);
	int initAtomsSize = atoms->getSize();
	atoms->resize(2);
	double dHalfCell = 0.5 * length / n;
	vector<double> transformV = { dHalfCell, dHalfCell, dHalfCell };
	for (int i = 0; i < initAtomsSize; i++) {
		vector<double> p = atoms->getPos(i);
		for (int k = 0; k < 3; k++) {
			p[k] += transformV[k];
		}
		atoms->setPos(i + initAtomsSize, p);
	}
	atoms->setCellLength(2.0 * dHalfCell);
	atoms->repeat(n);
	cout << "Chose BCC!" << endl;
}

void CellBuilder::buildFCCUnitCell(Atoms* atoms, double N, double length) {
	int n = static_cast<int>(N);
	int initAtomsSize = atoms->getSize();
	atoms->resize(4);
	double dHalfCell = 0.5 * length / n;
	vector<vector<double>> transformVs = {
		{ 0.0, dHalfCell, dHalfCell },
		{ dHalfCell, 0.0, dHalfCell },
		{ dHalfCell, dHalfCell, 0.0 }
	};
	for (int i = 0; i < initAtomsSize; i++) {
		for (int j = 0; j < 3; j++) {
			vector<double> p = atoms->getPos(i);
			for (int k = 0; k < 3; k++) {
				p[k] += transformVs[j][k];
			}
			atoms->setPos((j + (__int64)1) * initAtomsSize + i, p);
		}
	}
	atoms->setCellLength(2.0 * dHalfCell);
	atoms->repeat(n);
	cout << "Chose FCC!" << endl;
}