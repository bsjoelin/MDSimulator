// MDsimulator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include "Atoms.h"
#include "CellBuilder.h"
#include "VelocityManager.h"
#include "Ensemble.h"
#include "Analysis.h"
#include "dataType.h"
using namespace std;

const double kB = 1.38065e-23;
const double AVOGADRO = 6.022045e+23;
const string INFILE = "parameters.inp";
const string OUTFILE = "sim.out";


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

	struct dataT dataContainer;
	GetParameters(&dataContainer);
	Atoms atoms(dataContainer.nAtoms);
	InitializeSetup(&atoms, &dataContainer);
	Ensemble* ens = new NVE(&atoms, &dataContainer);
	AnalysisTools::LinearRegressor reg = AnalysisTools::LinearRegressor();

	double U, K, t = 0;
	logger << "t\tU\tK\tH" << endl;
	U = ens->calculate() * dataContainer.epsK;
	K = atoms.getEnergy() * dataContainer.epsK;
	logger << t << "\t" << U << "\t" << K << "\t" << K + U << endl;
	reg.addPoint(t, K + U);

	for (int i = 1; i <= dataContainer.simSteps; i++)
	{
		ens->update();
		U = ens->calculate() * dataContainer.epsK;
		K = atoms.getEnergy() * dataContainer.epsK;

		t += dataContainer.dt_ps;
		logger << t << "\t" << U << "\t" << K << "\t" << K + U << endl;
		reg.addPoint(t, K + U);
	}

	logger.close();
	cout << "dt = " << dataContainer.dt_ps << endl;
	cout << "a = " << reg.getSlope() << endl;
	cout << "b = " << reg.getIntersect() << endl;
	return 0;
}




void GetParameters(dataT *d) {
	ifstream inputFile(INFILE);
	if (inputFile.is_open()) {
		//This should be sent to a parser to allow better control of input files
		inputFile >> d->nAtoms >> d->simSteps >> d->dt_ps >> d->T >>
			d->rho >> d->P >> d->epsK >> d->sigma >> d->r_co;

		d->tau_s = 0.0;

		// Calculate values
		double mu = (Atoms::mass * Atoms::mass) / (2 * Atoms::mass) 
			/ (AVOGADRO * 1000.0);

		d->dt_s = pow(d->epsK * kB / (mu * pow(d->sigma, 2)), 0.5) 
			* 1.0e-12 * 1.0e+10 * d->dt_ps;

		// The reduced temperature is given by T_s = T * kB / eps = T / epsK
		d->T_s = d->T / d->epsK;

		// The reduced relaxation can be determined from the factor for dt_s
		d->tau_s_s = d->tau_s * d->dt_s / d->dt_ps;

		// Set integrator and potential
		d->IT = InteType::VELVERLET;
		d->PT = PotType::LJ;
		
		// Close the input file.
		inputFile.close();
	}
	else {
		throw invalid_argument("Data file unavailable!");
	}
}

void InitializeSetup(Atoms* a, dataT* d) {
	double rhoN = AVOGADRO / a->mass * d->rho * 1e-24 * pow(d->sigma, 3);
	CellBuilder::buildCell(a, rhoN, d->sigma);
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
