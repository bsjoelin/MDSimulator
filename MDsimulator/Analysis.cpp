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
