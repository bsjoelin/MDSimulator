#include "Parser.h"
#include "Ensemble.h"
#include "Potential.h"
#include "Integrator.h"

// Creates the InputParser with the alias matrix, and parses the input file.
// Then parses the values into the dataT object.
Parser::Parser(std::string f, dataT* d)
	: ip(aliasMatrix)
{
	ip.parseFile(f);
	parseData(d);
}

// empty destructor
Parser::~Parser() {}

// Parses all the necessary values into their dataT variable
void Parser::parseData(dataT* d) {
	int N = ip.getInt("N");
	if (N == -1) {
		std::cout << "Parameter 'N' (number of molecules) has to be given!"
			<< endl;
		exit(-1);
	}
	parseValue(&(d->nMolecules), "N");
	parseValue(&(d->apm), "apm");
	parseValue(&(d->simSteps), "steps");
	parseValue(&(d->mass), "mass");
	parseValue(&(d->T), "T");
	parseValue(&(d->rho), "rho");
	parseValue(&(d->P), "P");
	parseValue(&(d->dt_ps), "dt");
	parseValue(&(d->epsK), "epsK");
	parseValue(&(d->sigma), "sigma");
	parseValue(&(d->r_co), "r_c");
	parseValue(&(d->tau_s), "tau_s");
	parseValue(&(d->pos), "pos");
	parseValue(&(d->bonds), "bonds");
	parseValue(&(d->ks), "bks");
	parseValue(&(d->r_eqs), "r_eqs");
	parseValue(&(d->ET), "ens");
	parseValue(&(d->PT), "pot");
	parseValue(&(d->IT), "int");
}

// Following are all the functions for securely parsing a value into
// the dataT object, ensuring the type is correct, and if values
// were not given, the defaults take over.

void Parser::parseValue(int* vp, std::string key) {
	int val = ip.getInt(key);
	if (val != -1) {
		*vp = val;
	}
}

void Parser::parseValue(double* vp, std::string key) {
	double val = ip.getDouble(key);
	if (!(val + 1.0 < 0.0000001)) {
		*vp = val;
	}
}

void Parser::parseValue(std::string* vp, std::string key) {
	std::string val = ip.getString(key);
	if (val.compare("") != 0) {
		*vp = val;
	}
}

void Parser::parseValue(std::vector<double>* vp, std::string key) {
	*vp = ip.getVectorD(key);
}

void Parser::parseValue(std::vector<int>* vp, std::string key) {
	*vp = ip.getVectorI(key);
}

void Parser::parseValue(EnsType* vp, std::string key) {
	std::string val = ip.getString(key);
	if (val.compare("NVT") == 0 || val.compare("nvt") == 0) {
		*vp = EnsType::NVT;
	} else {
		*vp = EnsType::NVE;
	}
}

void Parser::parseValue(PotType* vp, std::string key) {
	/*std::string val = ip.getString(key);
	if (val.compare("LJ") == 0 || val.compare("lj") == 0) {
	*vp = PotType::LJ;
	}*/  // Future-proofing
	*vp = PotType::LJ;
}

void Parser::parseValue(InteType* vp, std::string key) {
	std::string val = ip.getString(key);
	if (val.compare("VELVERLET") == 0 || val.compare("velverlet") == 0
		|| val.compare("VelVerlet") == 0) {
		*vp = InteType::VELVERLET;
	} else {
		*vp = InteType::VERLET;
	}
}