#include "CellBuilder.h"
#include "math.h"

void CellBuilder::buildCell(Atoms* atoms, double density, double sigma) {
	int N = atoms->getSize();
	double n = pow(N, (1.0 / 3.0));
	double cellLength = pow(N / density, (1.0 / 3.0)) / sigma;
	double temp = 0.0;
	if (abs(modf(n + 0.0001, &temp)) < 0.001) {
		buildSCCell(temp, cellLength, atoms);
	}
	else if (abs(modf(n * bccFactor + 0.001, &temp)) < 0.01) {
		buildBCCCell(temp, cellLength, atoms);
	}
	else if (abs(modf(n * fccFactor + 0.001, &temp)) < 0.01) {
		buildFCCCell(temp, cellLength, atoms);
	}
	else {
		N = int(pow(floor(n), 3.0));
		cout << "Number of atoms is not a magic number! using "
			<< N << " instead." << endl;
		*atoms = Atoms(N);
		cellLength = pow(N / density, (1.0 / 3.0)) / sigma;
		buildSCCell(int(floor(n)), cellLength, atoms);
	}
	atoms->setCellLength(cellLength);
}

void CellBuilder::buildSCCell(double N, double length, Atoms * atoms) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				vector<double> r = { i * dCell, j * dCell, k * dCell};
				atoms->setPos(index++, r);
			}
		}
	}
	atoms->center();
	cout << "Chose SC!" << endl;
}

void CellBuilder::buildBCCCell(double N, double length, Atoms* atoms) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				vector<double> r = { i * dCell, j * dCell, k * dCell };
				atoms->setPos(index++, r);
				r = { (i + 0.5) * dCell, (j + 0.5) * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
			}
		}
	}
	atoms->center();
	cout << "Chose BCC!" << endl;
}

void CellBuilder::buildFCCCell(double N, double length, Atoms* atoms) {
	int n = static_cast<int>(N);
	double dCell = length / n;
	int index = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				vector<double> r = { i * dCell, j * dCell, k * dCell };
				atoms->setPos(index++, r);
				r = { (i + 0.5) * dCell, (j + 0.5) * dCell, k * dCell };
				atoms->setPos(index++, r);
				r = { i * dCell, (j + 0.5) * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
				r = { (i + 0.5) * dCell, j * dCell, (k + 0.5) * dCell };
				atoms->setPos(index++, r);
			}
		}
	}
	atoms->center();
	cout << "Chose FCC!" << endl;
}