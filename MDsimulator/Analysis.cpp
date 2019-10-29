#include "Analysis.h"
#include <iostream>

AnalysisTools::LinearRegressor::LinearRegressor()
{
}

void AnalysisTools::LinearRegressor::addPoint(double x, double y) {
	sumX += x;
	sumXX += x * x;
	sumY += y;
	sumXY += x * y;
	elements++;
}

double AnalysisTools::LinearRegressor::getSlope() {
	if (elements == 0) {
		return 0;
	}
	
	double averageX = sumX / elements;
	double averageY = sumY / elements;

	return (sumXY - elements * averageX * averageY) /
		(sumXX - elements * averageX * averageX);
}

double AnalysisTools::LinearRegressor::getIntersect() {
	if (elements == 0) {
		return 0;
	}

	double averageX = sumX / elements;
	double averageY = sumY / elements;

	return averageY - getSlope() * averageX;

}

void AnalysisTools::LinearRegressor::printSums() {
	std::cout << sumX << std::endl;
	std::cout << sumXX << std::endl;
	std::cout << sumXY << std::endl;
	std::cout << sumY << std::endl;
}
