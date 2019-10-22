
#ifndef _cellbuilder_h
#define _cellbuilder_h

#include "Atoms.h"
#include <iostream>
#include <cmath>
#include "math.h"

class CellBuilder
{
public:

	static void buildCell(Atoms* atoms, double density, double sigma);

private:
	
	static constexpr double bccFactor = 0.7937005259841;
	static constexpr double fccFactor = 0.6299605249474;

	static void buildSCCell(double N, double length, Atoms* atoms);
	static void buildBCCCell(double N, double length, Atoms* atoms);
	static void buildFCCCell(double N, double length, Atoms* atoms);
};

#endif // !_cellbuilder_h
