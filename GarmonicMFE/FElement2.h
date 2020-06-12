#pragma once
#include <vector>
#include <ostream>

struct Point2 {
	Point2() = default;
	~Point2() = default;

	Point2(const Point2& point) = default;
	Point2(Point2&& point) = default;

	Point2(double x, double y) : x(x), y(y) {};
	Point2(const std::pair<double, double>& copy) : x(copy.first), y(copy.second) {};

	Point2 operator=(const Point2& copy) {
		x = copy.x;
		y = copy.y;
		return *this;
	}

	double x;
	double y;
};

// Two-direction quadrange elem
class FElement2 {
public:
	FElement2() = default;
	~FElement2() = default;

	Point2 getPoint(size_t localIndex) const;
	void setPoint(const Point2& curPoint, size_t index);

	double getSigma() const;
	void setSigma(double value);

	//double getLambda() const;
	//void setLambda(double value);
	friend std::ostream& operator<<(std::ostream& out, const FElement2& grid);
private:

	Point2 vertic[4];
	
	//double lambda;
	double sigma;

	// Базисные функции сюда же в виде лямбд
};

