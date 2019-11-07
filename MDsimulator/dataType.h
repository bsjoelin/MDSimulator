#ifndef _datatype_h
#define _datatype_h

// Structure class to contain the parameters of the MD simulation
struct dataT {
	int nAtoms;			// Number of atoms
	int simSteps;		// Number of MD steps
	double T;			// Temperature [Kelvin]
	double rho;			// Density [g/cm^3]
	double P;			// Pressure [Pa]
	double dt_ps;		// delta time [ps]
	double epsK;		// epsilon/k_B [Kelvin] 
	double sigma;		// sigma [Angstrom]
	double r_co;		// Potential cut_off [Angstrom]
	double tau_s;		// Relaxation time for heat bath [ps]

	// Derived values
	double eps;			// epsilon [eV]

	// Dimensionless values
	double dt_s;		// Reduced timestep
	double T_s;			// Reduced temperature
	double tau_s_s;		// Reduced heat bath relaxation time 

	// Simulation type
	InteType IT;		// The integration scheme employed
	PotType PT;			// The potential employed
};

#endif // !_datatype_h

