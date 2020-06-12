#include "Gauss2.h"

const double Gauss2::csi[4] = { -1. / sqrt(3), -1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3) };
const double Gauss2::etta[4] = { -1. / sqrt(3), 1. / sqrt(3), -1. / sqrt(3), 1. / sqrt(3) };

Gauss2::Gauss2() {
	x1 = x2 = y1 = y2 = 0;
	jacobian = 0;
}

void Gauss2::setArea(double x1, double x2, double y1, double y2) {
	this->x1 = x1;
	this->x2 = x2;
	this->y1 = y1;
	this->y2 = y2;
	initArea();
}

// Integrate f, x=x1..x2, y=y1..y2, by Gauss2
double Gauss2::integrate(const std::function<double(double, double)>& f) const {
	double result = 0;
	for (size_t i = 0; i < 4; i++) {
		result += jacobian * f(convertKsiToX(csi[i]), convertEttaToY(etta[i]));
	}
	return result;
}

// Integrate mult(f1 * f2), x=x1..x2, y=y1..y2, by Gauss2
double Gauss2::integrate(const std::function<double(double, double)>& f1, const std::function<double(double, double)>& f2) const {
	double result = 0;
	for (size_t i = 0; i < 4; i++) {
		result += jacobian * f1(convertKsiToX(csi[i]), convertEttaToY(etta[i])) * f2(convertKsiToX(csi[i]), convertEttaToY(etta[i]));
	}
	return result;
}

bool Gauss2::isInit() const {
	return x1 != x2 && y1 != y2;
}

void Gauss2::initArea() {
	jacobian = (x2 - x1) * (y2 - y1) / 4;
}

double Gauss2::convertKsiToX(double ksi) const {
	return (x2 - x1) * (ksi + 1) / 2 + x1;
}
	
double Gauss2::convertEttaToY(double etta) const {
	return (y2 - y1) * (etta + 1) / 2 + y1;
}