// MDsimulator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include "Atoms.h"
#include "CellBuilder.h"
#include "VelocityManager.h"
#include "Ensemble.h"
using namespace std;

const double kB = 1.38065e-23;
const double AVOGADRO = 6.022045e+23;
const string INFILE = "parameters.inp";
const string OUTFILE = "sim.out";

struct dataT {
	int nAtoms, simSteps;  //number of atoms, number of MD steps.
	double T, rho, P, dt_ps;  //temperature, density, pressure.
	double epsK, sigma, r_co;  // epsilon/K, sigma, potential cut_off.
	double dt_s, T_s;  // Reduced timestep, reduced temperature
} dataContainer;

void GetParameters(dataT* data);
void InitializeSetup(Atoms* atoms, dataT* data);


int main()
{
	ofstream logger;
	logger.open(OUTFILE);
	if (!logger.is_open()) {
		cout << "Couldn't open output file. Exiting." << endl;
		return 0;
	}
	GetParameters(&dataContainer);
	Atoms atoms(dataContainer.nAtoms);
	InitializeSetup(&atoms, &dataContainer);
	Ensemble* ens = new NVE(&atoms, PotType::LJ, InteType::VERLET, dataContainer.dt_ps);
	double U, K;
	logger << "U\tK\tH" << endl;
	for (int i = 0; i < dataContainer.simSteps; i++)
	{
		U = ens->calculate();
		ens->update();

		K = atoms.getEnergy();
		logger << U << "\t" << K << "\t" << U + K << endl;
	}
	return 0;
}




void GetParameters(dataT *d) {
	ifstream inputFile(INFILE);
	if (inputFile.is_open()) {
		//This should be sent to a parser to allow better control of input files
		inputFile >> d->nAtoms >> d->simSteps >> d->dt_ps >> d->T >>
			d->rho >> d->P >> d->epsK >> d->sigma >> d->r_co;

		// do something
		double mu = (Atoms::mass * Atoms::mass) / (2 * Atoms::mass) 
			/ (AVOGADRO * 1000.0);

		d->dt_s = pow(d->epsK * kB / (mu * pow(d->sigma, 2)), 0.5) 
			* 1.0e-12 * 1.0e+10 * d->dt_ps;

		d->T_s = d->T / d->epsK;
	}
	else {
		throw invalid_argument("Data file unavailable!");
	}
	inputFile.close();
}

void InitializeSetup(Atoms* a, dataT* d) {
	CellBuilder::buildCell(a, d->rho, d->sigma);
	VelocityManager vm;
	vm.initializeVelocities(a, d->T_s);
}


// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
