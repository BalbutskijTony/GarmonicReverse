#include "StabilizedBiorderedGradient.h"

void StabilizedBiorderedGradient::setEndType(TypeEnd type) {

	switch (type)
	{
	case StabilizedBiorderedGradient::TypeEnd::BY_OMEGA:
		endCondition = [&]() {return abs(w) < criticalW; };
		break;
	case StabilizedBiorderedGradient::TypeEnd::BY_ERROR:
		endCondition = [&]() {return norm(r) / norm(*right) < criticalE; };
		break;
	case StabilizedBiorderedGradient::TypeEnd::ONLY_BY_ITER:
		endCondition = []() {return false; };
		break;
	default:
		break;
	}

}

void StabilizedBiorderedGradient::initZero(std::vector<double>& vect, size_t count) {
	if (!vect.empty()) vect.clear();
	for (size_t i = 0; i < count; i++)
		vect.push_back(0.0);
}

void StabilizedBiorderedGradient::initAlocate(size_t dimention) {
	initZero(r, dimention);
	initZero(rTild, dimention);
	initZero(p, dimention);
	initZero(v, dimention);
	initZero(t, dimention);
	initZero(s, dimention);

	if (x.empty()) initZero(x, dimention);

	roCur = 0;
	roLast = 0;
	a = 0;
	w = 0;
}

const std::vector<double>& StabilizedBiorderedGradient::solve(const std::vector<std::vector<double>>& matrix, const std::vector<double>& right) {
	(this->right) = &right;
	initAlocate(right.size());
	zeroIter(matrix, right);
	while (!isEnd()) {
		iterate(matrix, right);
	}
	return x;
}

void StabilizedBiorderedGradient::iterate(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	// 1.
	roCur = scal(rTild, r); 	
	// 2.
	b = roCur / roLast * a / w; 
	// 3.
	for (size_t i = 0; i < vec.size(); i++) {
		p[i] = r[i] + b * (p[i] - w * v[i]);
	}
	// 4.
	mult(matrix, p, v);
	// 5.
	double tmp;
	a = roCur / ((tmp = scal(rTild, v)) == 0 ? 1 : tmp); // странное поведение
	// 6.
	for (size_t i = 0; i < vec.size(); i++) {
		s[i] = r[i] - a * v[i];
	}
	// 7.
	mult(matrix, s, t);
	// 8.
	w = scal(t, s) / scal(t, t);

	for (size_t i = 0; i < vec.size(); i++) {
		// 9.
		x[i] = x[i] + w * s[i] + a * p[i];
		// 10.
		r[i] = s[i] - w * t[i];
	}
	roLast = roCur;
	curIter++;
}

bool StabilizedBiorderedGradient::isEnd() {
	return (curIter >= maxIter) || endCondition();
}

void StabilizedBiorderedGradient::zeroIter(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	mult(matrix, x, r);
	for (size_t i = 0; i < vec.size(); i++) {
		r[i] = vec[i] - r[i];
		rTild[i] = r[i];
	}
	roCur = roLast = a = w = 1;
	curIter = 0;
	if (0 == maxIter) maxIter = 1000;

}

//----------------------------------------------------------------------------------------------------------------------
void StabilizedBiorderedGradient::copy(std::vector<double>& to, const std::vector<double>& from) {
	for (size_t i = 0; i < from.size(); i++) {
		to[i] = from[i];
	}
}

void StabilizedBiorderedGradient::mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result) {
	for (size_t row = 0; row < vec.size(); row++) {
		result[row] = 0;
		for (size_t col = 0; col < vec.size(); col++) {
			result[row] += matrix[row][col] * vec[col];
		}
	}
}

std::vector<double> StabilizedBiorderedGradient::mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	std::vector<double> result;
	for (size_t row = 0; row < vec.size(); row++) {
		result.push_back(0.0);
		for (size_t col = 0; col < vec.size(); col++) {
			result[row] += matrix[row][col] * vec[col];
		}
	}
	return result;
}

double StabilizedBiorderedGradient::scal(const std::vector<double>& v1, const std::vector<double>& v2) const {
	double res = 0.0;
	for (size_t ind = 0; ind < v1.size(); ind++)
		res += v1[ind] * v2[ind];
	return res;
}

double StabilizedBiorderedGradient::norm(const std::vector<double>& v) {
	double result = 0;
	for (size_t i = 0; i < v.size(); i++)
		result += v[i] * v[i];
	return sqrt(result);
}

bool StabilizedBiorderedGradient::isCorrect() const {
	return endCondition();
}

void StabilizedBiorderedGradient::setMaxIter(size_t iterCount) {
	maxIter = iterCount;
}

void StabilizedBiorderedGradient::setCritW(double wCrit) {
	criticalW = wCrit;
}

void StabilizedBiorderedGradient::setCritRNorm(double norm) {
	criticalE = norm;
}

void StabilizedBiorderedGradient::setStart(const std::vector<double>& start) {
	if (!x.empty()) x.clear();
	
	for (size_t i = 0; i < start.size(); i++)
		x.push_back(start[i]);

}

size_t StabilizedBiorderedGradient::getIterCount() const {
	return curIter;
}