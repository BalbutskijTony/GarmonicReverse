#pragma once
#include <string>
#include "FEGrid2.h"

#include "LOS.h"
#include "StabilizedBiorderedGradient.h"

#include "Gauss2.h"
#include "BuilderGlobal.h"

class GarmonicTest
{
public:
	GarmonicTest();
	void setRealFunction(std::function<double(double, double)>& uSin, std::function<double(double, double)>& uCos);
	void setRightFunction(std::function<double(double, double)>& fSin, std::function<double(double, double)>& fCos);

	void createGrid(const std::string& gridFile, const std::string& xKoeff, const std::string& yKoeff);
	void setLogOut(std::ostream& out);
	void testing();
private:
	// Functions region START
	std::function<double(double, double)> uSin;
	std::function<double(double, double)> uCos;

	std::function<double(double, double)> fSin;
	std::function<double(double, double)> fCos;
	// Functions region END
	std::string name;
	Grid2 mainGrid;
	Grid2 fullGrid;
	double w, lambda, hi, sigma;
	FEGrid2 feGrid;
	BuilderGlobal<DIM_LOCAL> builder;
	StabilizedBiorderedGradient solver;

	std::ostream* logOut;
};

