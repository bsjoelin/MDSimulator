#ifndef _Atoms_h
#define _Atoms_h

#include <vector>

using namespace std;


class Atoms {
public:
	Atoms(int natoms);
	~Atoms();

	int getSize();
	vector<double> getPos(int i);
	vector<double> getVel(int i);
	void resize(int n);

	void setPos(int i, vector<double> r);
	void setVel(int i, vector<double> r);
	void setCellLength(double length);

	void center();
	void centerVel();
	void print();
	void printVel();

	static constexpr double mass = 39.948;

private:
	int nAtoms;
	double cellLength;
	vector<vector<double> > pos, oldPos, vel, oldVel;
};

#endif