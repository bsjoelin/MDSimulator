#ifndef _analysistools_h
#define _analysistools_h

class AnalysisTools
{
public:
	class LinearRegressor
	{
	public:
		LinearRegressor();
		
		void addPoint(double x, double y);
		double getSlope();
		double getIntersect();
		void printSums();

	private:
		double sumX = 0;
		double sumXX = 0;
		double sumXY = 0;
		double sumY = 0;
		int elements = 0;
	};

};

#endif // !_analysistools_h
