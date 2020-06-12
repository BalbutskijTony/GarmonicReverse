#pragma once
#include <functional>

class BasicFuncFabric
{
public:
	BasicFuncFabric();
	~BasicFuncFabric() = default;

	void setArea(double x1, double x2, double y1, double y2);

	std::function<double(double, double)> getBasic(size_t index) const;

private:
	double x1, x2, y1, y2;
	double hx, hy;

	std::function<double(double)> funcX1;
	std::function<double(double)> funcX2;
	std::function<double(double)> funcY1;
	std::function<double(double)> funcY2;

};

