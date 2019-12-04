#define _USE_MATH_DEFINES
#include "Analysis.h"
#include <iostream>

// Empty constructor, since all members are already initialized to zero
AnalysisTools::LinearRegressor::LinearRegressor(){}

// Adding a point adds to all the members and raises the element count by one
void AnalysisTools::LinearRegressor::addPoint(double x, double y) {
	sumX += x;
	sumXX += x * x;
	sumY += y;
	sumXY += x * y;
	elements++;
}

// Getting the slope uses a nice formula for calculating the slope from the sums
double AnalysisTools::LinearRegressor::getSlope() {
	// We can't have a slope, if there are no elements
	if (elements == 0) {
		return 0;
	}
	
	// Calculate averages
	double averageX = sumX / elements;
	double averageY = sumY / elements;

	// Simply return the result of the formula
	return (sumXY - elements * averageX * averageY) /
		(sumXX - elements * averageX * averageX);
}

// Getting the intersect uses  y = ax + b for calculating the intersect from
// the sums
double AnalysisTools::LinearRegressor::getIntersect() {
	// We can't have an intersect, if there are no elements
	if (elements == 0) {
		return 0;
	}

	// Calculate averages
	double averageX = sumX / elements;
	double averageY = sumY / elements;

	// y = ax + b => b = y - ax, but using averages
	return averageY - getSlope() * averageX;
}

// Simply print out the internal members to the console
void AnalysisTools::LinearRegressor::printSums() {
	std::cout << sumX << std::endl;
	std::cout << sumXX << std::endl;
	std::cout << sumXY << std::endl;
	std::cout << sumY << std::endl;
}


AnalysisTools::RadDistribFunc::RadDistribFunc(Atoms* a, dataT* d)
	: hist(0, 0)
{
	rMax = a->getCellLength() / 2.0;
	int nbins = static_cast<int>(rMax / 0.02);
	dr = rMax / nbins;
	hist.resize(nbins);
	rhoN = d->rhoN;
	atoms = a;
}

AnalysisTools::RadDistribFunc::~RadDistribFunc() {
	vector<int>().swap(hist);		// Release memory from hist
}

void AnalysisTools::RadDistribFunc::update() {
	vector<vector<double>> r = atoms->getDistances();
	for (int i = 0; i < atoms->getSize()-1; i++) {
		for (int j = i+1; j < atoms->getSize(); j++) {
			if (r[i][j] > rMax) {
				continue;
			}
			int index = static_cast<int>( r[i][j] / dr);
			hist[index] += 2;
		}
	}
	nt++;
	vector<vector<double>>().swap(r);	// Release memory from r
}

vector<vector<double>> AnalysisTools::RadDistribFunc::getRDF() {
	double prefactor = atoms->getSize() * static_cast<double>(nt) 
		* 4.0 * M_PI * rhoN * pow(dr, 3.0) / 3.0;
	vector<vector<double>> RDF(hist.size(), vector<double>(2, 0));
	for (int i = 0; i < hist.size(); i++) {
		RDF[i][0] = (i + 0.5) * dr;
		RDF[i][1] = hist[i] / prefactor / (pow(i + 1, 3.0) - pow(i, 3.0));
	}
	return RDF;
}


AnalysisTools::Diffusion::Diffusion(Atoms* a, dataT* d)
	: op(0, vector<double>(0, 0))
{
	atoms = a;
	data = d;
}

AnalysisTools::Diffusion::~Diffusion() {

}

void AnalysisTools::Diffusion::start(double time) 
{
	for (int i = 0; i < atoms->getSize(); i++)
	{
		op.push_back(atoms->getPos(i));
	}
	startt = time;
}

double AnalysisTools::Diffusion::getDiffu(double time)
{
	if (op.size() == 0)
	{
		return 0;
	}
	for (int i = 0; i < atoms->getSize(); i++)
	{
		vector<double> np = atoms->getPos(i);
		double r =  0.0;
		for (int k = 0; k < 3; k++)
		{
			r += pow(np[k] - op[i][k], 2.0);
		}
		msd += r;
	}
	endt = (time - startt);
	return msd / (6.0 * endt * atoms->getSize());
}
