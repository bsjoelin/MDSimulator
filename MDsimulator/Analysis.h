#ifndef _analysistools_h
#define _analysistools_h
#include <vector>
#include "Atoms.h"
#include "dataType.h"


// Container class for analysis function/tools
class AnalysisTools
{
public:
	// A linear regression class
	class LinearRegressor
	{
	public:
		// Constructor
		LinearRegressor();
		
		// Add a (x, y) point for the regressor
		void addPoint(double x, double y);
		// Get the slope of the linear regression. Can be called at any point
		double getSlope();
		// Get the intersect of the linear regression. Can be called at any point
		double getIntersect();
		// Print the private members to the console
		void printSums();

	private:
		double sumX = 0;	// sum of x
		double sumXX = 0;	// sum of x * x
		double sumXY = 0;	// sum of x * y
		double sumY = 0;	// sum of y
		int elements = 0;	// the number of points in the regression
	};


	// Radial Distribution Class
	class RadDistribFunc
	{
	public:
		// Constructor
		RadDistribFunc(Atoms* atoms, dataT* data);
		// Destructor
		~RadDistribFunc();
		// Build histogram
		void update();
		// Get radial distribution function
		vector<vector<double>> getRDF();

	private:
		double dr;			// Delta r is the size of a bin
		int nt = 0;				// Counter for number of recorded configs. 
		vector<int> hist;	// Histogram
		double rhoN;		// Number density
		Atoms* atoms;
		double rMax;
	};
};

#endif // !_analysistools_h
