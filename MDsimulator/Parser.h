#ifndef _parser_h
#define _parser_h

#include "InputParser.h"
#include "dataType.h"
#include <iostream>

// Forward declarations
enum class EnsType;
enum class PotType;
enum class InteType;

class Parser
{
public:
	// Constructor loads the input file and loads the values into the dataT struct
	Parser(std::string filename, dataT* dataContainer);
	virtual ~Parser();

private:
	void parseData(dataT* dataContainer);

	// Overrides for parse a value of a specific type into the dataT struct.
	// All have in-built error handling (some inplicit).
	void parseValue(int* valptr, std::string key);
	void parseValue(double* valptr, std::string key);
	void parseValue(std::string* valptr, std::string key);
	void parseValue(std::vector<double>* valptr, std::string key);
	void parseValue(std::vector<int>* valptr, std::string key);
	void parseValue(EnsType* valptr, std::string key);
	void parseValue(PotType* valptr, std::string key);
	void parseValue(InteType* valptr, std::string key);

	// All the keywords with their associated aliases
	std::vector<std::vector<std::string>> aliasMatrix{
		{"N", "nmols"},
		{"apm", "atoms_per_mol"},
		{"steps", "simSteps"},
		{"mass"},
		{"dt", "timestep"},
		{"T", "temperature"},
		{"P", "pressure"},
		{"rho", "density"},
		{"epsK", "epsilon_in_K"},
		{"sigma"},
		{"r_c", "cutoff"},
		{"tau_s", "relaxation_time"},
		{"ens", "Ensemble"},
		{"pot", "Potential"},
		{"int", "Integrator"},
		{"pos", "positions"},
		{"bonds", "bond_pairs"},
		{"bks", "bond_constants"},
		{"r_eqs", "bond_eq_distances"}
	};

	// An Input Parser to parse the input through
	InputParser ip; // has to be at the bottom.
};

#endif // !_parser_h
