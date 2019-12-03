#include "Parser.h"
#include "Ensemble.h"
#include "Potential.h"
#include "Integrator.h"

Parser::Parser(std::string f, dataT* d)
	: ip(aliasMatrix)
{
	ip.parseFile(f);
	parseData(d);
}

// empty destructor
Parser::~Parser() {}


void Parser::parseData(dataT* d) {
	int N = ip.getInt("N");
	if (N == -1) {
		throw std::invalid_argument(
			"Parameter 'N' (number of molecules) has to be given!");
	}
	parseValue(&(d->apm), "apm");
	d->nAtoms = d->apm * N;
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