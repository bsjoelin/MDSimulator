// MDsimulator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>

#include <iomanip>

#include "Atoms.h"
#include "CellBuilder.h"
#include "VelocityManager.h"
#include "Ensemble.h"
#include "Analysis.h"
#include "dataType.h"
using namespace std;

// Define important constants
const double kB = 1.38064852e-23;  // Boltzmann's constant
const double kBeV = 8.6173303360e-5;  //Boltzmann's constant in eV/K
const double AVOGADRO = 6.022045e+23;  // Avogadro's constant
const string INFILE = "parameters.inp";  // Name of the input file - should be sysarg at some point.
const string OUTFILE = "sim.out";  // Name of output file - should be sysarg at some point.

// Function prototypes for main
void GetParameters(dataT* data);
void InitializeSetup(Atoms* atoms, dataT* data);
void saveXYZ(Atoms* atoms, dataT* data, string out);


// Main program execution routine
int main()
{
	// Create logger/output file as a stream and open it
	ofstream logger;
	logger.open(OUTFILE);
	
	// Make sure that the file was opened properly
	if (!logger.is_open()) {
		cout << "Couldn't open output file. Exiting." << endl;
		return 0;  // End the program, if the logger wasn't opened
	}

	// Create the data container and populate it
	dataT dataContainer;
	GetParameters(&dataContainer);
	// Initialize the Atoms object
	Atoms atoms(dataContainer.nAtoms);
	// Create the setup of the initial system
	InitializeSetup(&atoms, &dataContainer);

	// Create an Ensemble object, passing the Atoms object and data container
	Ensemble* ens = new NVE(&atoms, &dataContainer); // should take the Ensemble type as parameter
	
	// Initialize a linear regressor to take care of calculating the deviation in
	// the (extended) Hamiltonian
	AnalysisTools::LinearRegressor reg = AnalysisTools::LinearRegressor();

	// Initialize the potential and kinetic energy, and the time
	double U, K, Hx, t = 0;

	// Log the header
	logger << "t\tU\tK\tHx\tH" << endl;
	// Calculate and log the initial values
	U = ens->calculate() * dataContainer.eps;  // in eV
	K = atoms.getEnergy() * dataContainer.eps;  // in eV
	logger << t << "\t" << U << "\t" << K << "\t" << 0.0 << "\t" << K + U << endl;

	// Add the first point to the regressor
	reg.addPoint(t, K + U);  // Hx = 0 in the start

	// The MD loop of the program
	for (int i = 1; i <= dataContainer.simSteps; i++)
	{
		// Make the Ensemble update the positions and velocities of the Atoms object
		Hx = ens->update() * dataContainer.eps;
		// Calculate energies
		U = ens->calculate() * dataContainer.eps;
		K = atoms.getEnergy() * dataContainer.eps;

		// Log the time and energies
		t += dataContainer.dt_ps;  // actual time
		logger << t << "\t" << U << "\t" << K << "\t" << Hx << "\t" 
			<< K + U + Hx << endl;
		
		// Add the time-energy point to the regressor
		reg.addPoint(t, K + U + Hx);
	}

	// Close the logger and write regression data to the console
	logger.close();
	cout << "dt = " << dataContainer.dt_ps << endl;
	cout << "a = " << reg.getSlope() << endl;
	cout << "b = " << reg.getIntersect() << endl;
	cout << "P = " << ens->getPressure() * dataContainer.epsK * kB
		/ pow(dataContainer.sigma, 3.0) * 1e30 << " Pa" << endl;

	return 0;  // End program execution
}


// Function for population the data container from the input file
void GetParameters(dataT *d) {
	// Create input stream
	ifstream inputFile(INFILE);
	if (inputFile.is_open()) {
		// Stream the values of the input file into the correct values
		// (This should be sent to a parser to allow better control of input files)
		inputFile >> d->nAtoms >> d->simSteps >> d->dt_ps >> d->T >>
			d->rho >> d->P >> d->epsK >> d->sigma >> d->r_co >> d->tau_s;

		// Calculate reduced parameters
		double mu = (Atoms::mass * Atoms::mass) / (2 * Atoms::mass) 
			/ (AVOGADRO * 1000.0);

		d->dt_s = pow(d->epsK * kB / (mu * pow(d->sigma, 2)), 0.5) 
			* 1.0e-12 * 1.0e+10 * d->dt_ps;

		// The reduced temperature is given by T_s = T * kB / eps = T / epsK
		d->T_s = d->T / d->epsK;

		d->eps = d->epsK * kBeV;

		// The reduced relaxation can be determined from the factor for dt_s
		d->tau_s_s = d->tau_s * d->dt_s / d->dt_ps;

		// Set integrator and potential
		d->IT = InteType::VERLET;
		d->PT = PotType::LJ;
		
		// Close the input file.
		inputFile.close();
	}
	else {
		// End the program in error, if the inpu file is unavailable
		throw invalid_argument("Data file unavailable!");
	}
}

// Function for creating the initial configuration of the system
void InitializeSetup(Atoms* a, dataT* d) {
	// Calculate number density in m^{-3}
	double rhoN = AVOGADRO / a->mass * d->rho * 1e-24 * pow(d->sigma, 3);
	// Build the cell (with a static call)
	CellBuilder::buildCell(a, rhoN);
	// Initialize the velocities
	VelocityManager::initializeVelocities(a, d->T_s);
}


void saveXYZ(Atoms* a, dataT* d, string out) {
	ofstream outfile(out + ".xyz");
	if (outfile.is_open()) {
		outfile << a->getSize() << endl << endl;
		for (int i = 0; i < a->getSize(); i++) {
			outfile << setw(3) << "Ar";
			vector<double> q = a->getPos(i);
			for (int j = 0; j < 3; j++) {
				outfile << setw(15) << fixed << setprecision(5) 
					<< q[j] * d->sigma;
			}
			outfile << endl;
		}

		outfile.close();
	} else {
		throw invalid_argument("Couldn't open output file");
	}
}