#ifndef _datatype_h
#define _datatype_h

enum class InteType;
enum class PotType;
enum class EnsType;

// Structure class to contain the parameters of the MD simulation
struct dataT {
	int nAtoms = 1;			// Number of atoms
	int apm = 1;			// Atoms per molecules
	int simSteps = 1;		// Number of MD steps
	double T = 273.15;		// Temperature [Kelvin]
	double rho = 1.0;		// Density [g/cm^3]
	double mass = 1.0;		// Mass per atom [amu]
	double P = 101325;		// Pressure [Pa]
	double dt_ps = 0.001;	// delta time [ps]
	double epsK = 100;		// epsilon/k_B [Kelvin] 
	double sigma = 2.5;		// sigma [Angstrom]
	double r_co = 0.0;		// Potential cut_off [Angstrom]
	double tau_s = 0.0;		// Relaxation time for heat bath [ps]
	std::vector<double> pos{};		// The initial positions for the atoms
	std::vector<double> bonds{};	// The bonding pairs
	std::vector<double> ks{};		// The bonding force constants [eV/Angstrom^2]
	std::vector<double> r_eqs{};	// The equilibrium distances [Angstrom]

	// Derived values
	double eps = 0;			// epsilon [eV]

	// Dimensionless values
	double dt_s = 0;		// Reduced timestep
	double rhoN = 0;		// Reduced number density
	double T_s = 0;			// Reduced temperature
	double tau_s_s = 0;		// Reduced heat bath relaxation time 

	// Simulation type
	EnsType ET = EnsType(0);		// The ensemble type employed
	InteType IT = InteType(0);		// The integration scheme employed
	PotType PT = PotType(0);		// The potential employed
};

#endif // !_datatype_h

