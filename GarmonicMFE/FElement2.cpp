#include "FElement2.h"

Point2 FElement2::getPoint(size_t localIndex) const {
	return vertic[localIndex];
}

double FElement2::getSigma() const {
	return sigma;
}

void FElement2::setSigma(double value) {
	sigma = value;
}

//double FElement2::getGamma() const {
//	return gamma;
//}
//
//void FElement2::setGamma(double value) {
//	gamma = value;
//}
//
//double FElement2::getLambda() const {
//	return lambda;
//}
//
//void FElement2::setLambda(double value) {
//	lambda = value;
//}

void FElement2::setPoint(const Point2& curPoint, size_t index) {
	if (index > 3) return;

	vertic[index] = curPoint;
}

std::ostream& operator<<(std::ostream& out, const FElement2& grid) {
	for (size_t i = 0; i < 4; i++) {
		out << "P" << i + 1 << " = ( " << grid.vertic[i].x << ", " << grid.vertic[i].y << ")\t";
	}
	out << std::endl;
	return out;
}