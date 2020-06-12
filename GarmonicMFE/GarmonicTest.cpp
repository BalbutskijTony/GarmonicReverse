#include "GarmonicTest.h"

GarmonicTest::GarmonicTest(){
	logOut = nullptr;
}

void GarmonicTest::setRealFunction(std::function<double(double, double)>& uSin, std::function<double(double, double)>& uCos) {
	this->uSin = uSin;
	this->uCos = uCos;
}

void GarmonicTest::setRightFunction(std::function<double(double, double)>& fSin, std::function<double(double, double)>& fCos) {
	this->fSin = fSin;
	this->fCos = fCos;
}

void GarmonicTest::createGrid(const std::string& gridFile, const std::string& xKoeff, const std::string& yKoeff) {
	mainGrid.fromFile(gridFile);
	fullGrid = mainGrid.procreate(xKoeff, yKoeff);
	// ERROR!!!!
	//feGrid.constructFEGrid(fullGrid);

	if (nullptr != logOut) {
		*logOut << "MAIN Grid" << std::endl << mainGrid;
		*logOut << "FULL Grid" << std::endl << fullGrid;
		*logOut << "Finite Element Grid " << std::endl << feGrid;
	}
}

void GarmonicTest::setLogOut(std::ostream& out) {
	logOut = &out;
}

void GarmonicTest::testing() {
	builder.setFEGrid(feGrid);
	builder.setFullGrid(fullGrid);
	builder.setKernM(kernMLocal);
	builder.setKernG(kernGLocalXY, kernGLocalYX);

	builder.initGlobalMatrix();
	builder.initGlobalPAndC();

	builder.fillGlobalMatrix();
	// ERROR!!!
	//builder.fillGlobalPAndC(w, lambda, hi, sigma);

	//builder.setFSin(fSin);
	//builder.setFCos(fCos);

	// Not finished
}