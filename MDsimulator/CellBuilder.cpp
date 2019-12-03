#include "CellBuilder.h"
#include <iostream>
#include <cmath>
#include "math.h"

// The static buildCell() function determines the correct lattice sytem to
// build and sends the necessary parameters on to the proper functions
void CellBuilder::buildCell(Atoms* atoms, double density) {
	int N = atoms->getSize();  // the number of atoms
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
		buildSCCell(atoms, temp, cellLength);
	}
	else if (abs(modf(n * bccFactor + 0.001, &temp)) < 0.01) {
		buildBCCCell(atoms, temp, cellLength);
	}
	else if (abs(modf(n * fccFactor + 0.001, &temp)) < 0.01) {
		buildFCCCell(atoms, temp, cellLength);
	}
	// If c is neither 1, 2 or 4, we can assume, that the number is not a
	// magic number, and we just round N down to the nearest sc system.
	else {
		N = int(pow(floor(n), 3.0));
		cout << "Number of atoms is not a magic number! using "
			<< N << " instead." << endl;

		// Delete the old Atoms object and create a new with the proper
		// number of elements
		double mass = atoms->getMass();
		atoms->~Atoms();
		*atoms = Atoms(N, mass);
		cellLength = pow(N / density, (1.0 / 3.0));  // new cell size
		buildSCCell(atoms, n + 0.1, cellLength);
	}
	// Give the cell size to the Atoms object
	atoms->setCellLength(cellLength);
}

// All the following functions work by:
// 1) Determine the number of subcells (n^3)
// 2) Calculate the side length of each subcell
// 3) Place the atoms in each subcell as the lattice type prescribes

void CellBuilder::buildSCCell(Atoms * atoms, double N, double length) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)  // loop through x
	{
		for (int j = 0; j < n; j++)  // loop through y
		{
			for (int k = 0; k < n; k++)  // loop through z
			{
				// Place atom at (0, 0, 0) in each subcell
				vector<double> r = { i * dCell, j * dCell, k * dCell};
				atoms->setPos(index++, r);
			}
		}
	}
	// Center the atoms
	atoms->center();
	cout << "Chose SC!" << endl;
}

void CellBuilder::buildBCCCell(Atoms* atoms, double N, double length) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)  // loop through x
	{
		for (int j = 0; j < n; j++)  // loop through y
		{
			for (int k = 0; k < n; k++)  // loop through z
			{
				// Place atom at (0, 0, 0) in each subcell
				vector<double> r = { i * dCell, j * dCell, k * dCell };
				atoms->setPos(index++, r);
				// Place atom at (0.5, 0.5, 0.5) in each subcell (center)
				r = { (i + 0.5) * dCell, (j + 0.5) * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
			}
		}
	}
	// Center the atoms
	atoms->center();
	cout << "Chose BCC!" << endl;
}

void CellBuilder::buildFCCCell(Atoms* atoms, double N, double length) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)  // loop through x
	{
		for (int j = 0; j < n; j++)  // loop through y
		{
			for (int k = 0; k < n; k++)  // loop through z
			{
				// Place atom at (0, 0, 0) in each subcell
				vector<double> r = { i * dCell, j * dCell, k * dCell };
				atoms->setPos(index++, r);
				// Place atom at (0.5, 0.5, 0) in each subcell (xy-plane)
				r = { (i + 0.5) * dCell, (j + 0.5) * dCell, k * dCell };
				atoms->setPos(index++, r);
				// Place atom at (0.5, 0, 0.5) in each subcell (xz-plane)
				r = { i * dCell, (j + 0.5) * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
				// Place atom at (0, 0.5, 0.5) in each subcell (yz-plane)
				r = { (i + 0.5) * dCell, j * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
			}
		}
	}
	// Center the atoms
	atoms->center();
	cout << "Chose FCC!" << endl;
}