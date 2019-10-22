#ifndef _Atoms_h
#define _Atoms_h

#include <vector>
#include <iostream>

using namespace std;


class Atoms {
public:
	Atoms(int natoms);
	~Atoms();

	void center();
	void centerVel();
	void print(bool printAverage);
	void printVel();

	int getSize();
	vector<double> getPos(int i);
	vector<double> getVel(int i);
	double getEnergy();

	void setPos(int i, vector<double> r);
	void setVel(int i, vector<double> r);
	void setCellLength(double length);

	static constexpr double mass = 39.948;

private:
	int nAtoms;
	double cellLength;
	vector<vector<double> > pos, vel;
};

#endif