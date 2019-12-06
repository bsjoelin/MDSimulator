#ifndef _cellbuilder_h
#define _cellbuilder_h

#include "Atoms.h"

// Static class used for building the atomic lattice. Automatically detects the
// lattice type as either simple cubic, body centered cubic or face-centered cubic
// depending on the number of atoms. If the number fits none of these, the number 
// of atoms will be rounded down to fit the nearest simple cubic system.
class CellBuilder
{
public:
	// Static function for building the atomic system. Takes the Atoms object
	// to populate (as a pointer) and the number density of the system.
	static void buildCell(Atoms* atoms, int nMolecules, double density);

private:
	// Constants for the detection algorithm
	static constexpr double bccFactor = 0.7937005259841;  // cube root of 1/2
	static constexpr double fccFactor = 0.6299605249474;  // cube root of 1/4

	// Functions that build unit cell of the given system
	static void BuildSCUnitCell(Atoms* atoms, double N, double length);
	static void buildBCCUnitCell(Atoms* atoms, double N, double length);
	static void buildFCCUnitCell(Atoms* atoms, double N, double length);
};

#endif // !_cellbuilder_h
