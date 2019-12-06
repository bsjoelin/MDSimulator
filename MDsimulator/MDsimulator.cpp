// MDsimulator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Atoms.h"
#include "CellBuilder.h"
#include "VelocityManager.h"
#include "Ensemble.h"
#include "dataType.h"
#include "Analysis.h"
#include "Parser.h"
using namespace std;

// Define important constants
const double kB = 1.38064852e-23;  // Boltzmann's constant
const double kBeV = 8.6173303360e-5;  //Boltzmann's constant in eV/K
const double AVOGADRO = 6.022045e+23;  // Avogadro's constant
const string INFILE = "params.in";  // Name of the input file - should be sysarg at some point.
const string OUTFILE = "sim.out";  // Name of output file - should be sysarg at some point.

// Function prototypes for main
void GetParameters(dataT* data);
void InitializeSetup(Atoms* atoms, dataT* data);
void saveXYZ(Atoms* atoms, dataT* data, string out);


// Main program execution routine
int main()
{
	// Create logger/output file as a stream and open it
	ofstream logger(OUTFILE);
	
	// Make sure that the file was opened properly
	if (!logger.is_open()) {
		cout << "Couldn't open output file. Exiting." << endl;
		return -1;  // End the program, if the logger wasn't opened
	}

	// Create the data container and populate it
	dataT dataContainer;
	GetParameters(&dataContainer);

	// Initialize the Atoms object
	Atoms atoms(dataContainer.apm, dataContainer.mass);
	// Create the setup of the initial system
	InitializeSetup(&atoms, &dataContainer);
	
	// Create an Ensemble object, passing the Atoms object and data container
	Ensemble* ens = Ensemble::createEnsemble(&atoms, &dataContainer); // should take the Ensemble type as parameter

	// Initialize a linear regressor to take care of calculating the deviation in
	// the (extended) Hamiltonian
	AnalysisTools::LinearRegressor reg = AnalysisTools::LinearRegressor();
	AnalysisTools::RadDistribFunc rdf = AnalysisTools::RadDistribFunc(&atoms, &dataContainer);
	AnalysisTools::Diffusion dico = AnalysisTools::Diffusion(&atoms, &dataContainer);

	// Initialize the potential and kinetic energy, pressure and the time
	double U, K, Hx, t = 0;
	double avPressure = 0;

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
		if (i > 10000) {
			rdf.update();
			avPressure += ens->getPressure();
		}

		// Register start time and positions for self-diffusion
		if (i == 10001)	{
			dico.start(t * dataContainer.dt_s / dataContainer.dt_ps);
		}
	}

	// Calculate average pressure for the steps after 10000
	avPressure = avPressure/(dataContainer.simSteps - 10000.0);

	// Close the logger and write regression data to the console
	logger.close();
	cout << "dt = " << dataContainer.dt_ps << endl;
	cout << "a = " << reg.getSlope() << " eV/ps" << endl;
	cout << "b = " << reg.getIntersect() << " eV" << endl;
	cout << "p = " << ens->getPressure() * dataContainer.epsK * kB
		/ pow(dataContainer.sigma, 3.0) * 1e30 << " Pa" << endl;
	cout << "Z = " << avPressure/dataContainer.rhoN
		/dataContainer.T_s << endl;
	cout << "D = " << dico.getDiffu(t* dataContainer.dt_s 
		/ dataContainer.dt_ps) * pow(dataContainer.sigma, 2.0) * 1e-8
		<< " m^2/s" << endl;
	
	vector<vector<double>> graph = rdf.getRDF();
	ofstream rdfgraph("rdf.txt");
	
	// Make sure that the RDF file was opened properly
	if (!rdfgraph.is_open()) {
		cout << "Couldn't open rdf output file. Exiting." << endl;
		return -1;  // End the program, if the logger wasn't opened
	}
	// Print radial distribution function to rdf
	rdfgraph << "r" << "\t" << "g_r" << endl;
	for (vector<double> c : graph) {
		rdfgraph << c[0] << "\t" << c[1] << endl;
	}
	rdfgraph.close();

	saveXYZ(&atoms, &dataContainer, "fred");

	return 0;  // End program execution
}


// Function for population the data container from the input file
void GetParameters(dataT *d) {
	// Parse the input file
	Parser ps(INFILE, d);

	// Calculate reduced parameters
	double mu = (d->mass * d->mass) / (2 * d->mass) 
		/ (AVOGADRO * 1000.0);

	d->dt_s = pow(d->epsK * kB / (mu * pow(d->sigma, 2)), 0.5) 
		* 1.0e-12 * 1.0e+10 * d->dt_ps;

	// The reduced temperature is given by T_s = T * kB / eps = T / epsK
	d->T_s = d->T / d->epsK;

	d->eps = d->epsK * kBeV;

	// The reduced relaxation can be determined from the factor for dt_s
	d->tau_s_s = d->tau_s * d->dt_s / d->dt_ps;

	// Reduced bond distances and force constants
	double redFact = d->sigma * d->sigma / d->eps;
	for (double &k : d->ks) {
		k *= redFact;
	}

	for (double &r : d->r_eqs) {
		r /= d->sigma;
	}
}

// Function for creating the initial configuration of the system
void InitializeSetup(Atoms* a, dataT* d) {
	if (d->pos.size() != d->apm * 3.0) {
		string em = "Not enough positions were given! Found: "
			+ to_string(d->pos.size()) + " , but expected: "
			+ to_string(d->apm * (__int64)3);
		cout << em << endl;
		exit(-1);
	}

	for (int i = 0; i < d->apm; i++) {
		a->setPos(i, vector<double>{
			d->pos[(__int64)3 * i] / d->sigma,
			d->pos[(__int64)3 * i + (__int64)1] / d->sigma,
			d->pos[(__int64)3 * i + (__int64)2] / d->sigma
		});
	}

	if (d->bonds.size() % 2 == 1) {
		cout << "Bonds not given as pairs" << endl;
		exit(-1);
	}
	a->setBonds(d->bonds, d->ks, d->r_eqs);

	// Calculate number density in m^{-3}
	double rhoN = AVOGADRO / (d->apm * d->mass) * d->rho * 1e-24 * pow(d->sigma, 3);
	d->rhoN = rhoN;
	// Build the cell (with a static call)
	CellBuilder::buildCell(a, d->nMolecules, rhoN);
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
		cout << "Couldn't open output file" << endl;
		exit(-1);
	}
}