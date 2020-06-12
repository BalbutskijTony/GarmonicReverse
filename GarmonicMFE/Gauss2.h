#pragma once
#include <math.h>
#include <functional>

class Gauss2
{
public:
	Gauss2();
	~Gauss2() = default;

	void setArea(double x1, double x2, double y1, double y2);
	double integrate(const std::function<double(double, double)>& f) const;
	double integrate(const std::function<double(double, double)>& f1, const std::function<double(double, double)>& f2) const;

	bool isInit() const;
private:
	void initArea();
	double convertKsiToX(double ksi) const;
	double convertEttaToY(double etta) const;

	double x1, x2, y1, y2;
	double jacobian;

	static const double csi[4];
	static const double etta[4];
};
