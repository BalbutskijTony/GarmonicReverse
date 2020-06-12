#include "BasicFuncFabric.h"

BasicFuncFabric::BasicFuncFabric() {
	x1 = x2 = y1 = y2 = 0.0; 
	hx = hy = 0.0;
}

void BasicFuncFabric::setArea(double x1, double x2, double y1, double y2) {
	this->x1 = x1;
	this->x2 = x2;
	this->y1 = y1;
	this->y2 = y2;
	hx = x2 - x1;
	hy = y2 - y1;

	funcX1 = [=](double x) {return (x2 - x) / hx; };
	funcX2 = [=](double x) {return (x - x1) / hx; };
	funcY1 = [=](double y) {return (y2 - y) / hy; };
	funcY2 = [=](double y) {return (y - y1) / hy; };
}

// bi-linear index=0..3
std::function<double(double, double)> BasicFuncFabric::getBasic(size_t index) const {
	if (0 == index) return [=](double x, double y) {return funcX1(x) * funcY1(y); };
	if (1 == index) return [=](double x, double y) {return funcX2(x) * funcY1(y); };
	if (2 == index) return [=](double x, double y) {return funcX1(x) * funcY2(y); };
	if (3 == index) return [=](double x, double y) {return funcX2(x) * funcY2(y); };
}