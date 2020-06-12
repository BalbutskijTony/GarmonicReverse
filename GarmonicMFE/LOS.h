#pragma once
#include <math.h>
#include <ostream>
#include <optional>
#include <vector>

template<size_t precision, size_t maxIter> class LOS
{
public:
	LOS(/*std::optional<std::ostream&> out*/) {};

	std::vector<double> solve(const std::vector<std::vector<double>>& matrix, const std::vector<double>& right);
	void setStart(const std::vector<double>& start);
	bool isCorrect() const;
private:
	/*std::optional<std::ostream&> out;*/
	std::vector<double> startVect;

	// ------------------------------
	std::vector<double> z;
	std::vector<double> p;
	double a, b;
	std::vector<double> result;
	std::vector<double> r;
	// temp A*r_k
	std::vector<double> ar; 
	double curDelt;
	size_t curIter;
	// ------------------------------
	void copy(std::vector<double>& to, std::vector<double>& from);

	void mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result);
	std::vector<double> mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
	double scal(const std::vector<double>& v1, const std::vector<double>& v2) const;

	void initStart(size_t dimention);
	void initAlocate(size_t dimention);
	void iterate(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
	bool isEnd();
	void zeroIter(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
};

template<size_t precision, size_t maxIter>
std::vector<double> LOS<precision, maxIter>::solve(const std::vector<std::vector<double>>& matrix, const std::vector<double>& right) {
	initStart(right.size());
	initAlocate(right.size());
	zeroIter(matrix, right);
	while (!isEnd()) {
		iterate(matrix, right);
	}
	return result;
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::initAlocate(size_t dimention) {
	z.clear();
	p.clear();
	result.clear();
	r.clear();
	ar.clear();

	for (size_t i = 0; i < dimention; i++) {
		z.push_back(0.0);
		p.push_back(0.0);
		result.push_back(0.0);
		r.push_back(0.0);
		ar.push_back(0.0);
	}

	curIter = 0;
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result) {
	for (size_t row = 0; row < vec.size(); row++) {
		result[row] = 0;
		for (size_t col = 0; col < vec.size(); col++) {
			result[row] += matrix[row][col] * vec[col];
		}
	}
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::copy(std::vector<double>& to, std::vector<double>& from) {
	for (size_t i = 0; i < from.size(); i++) {
		to[i] = from[i];
	}
}

template<size_t precision, size_t maxIter>
std::vector<double> LOS<precision, maxIter>::mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	std::vector<double> result;
	for (size_t row = 0; row < vec.size(); row++) {
		result.push_back(0.0);
		for (size_t col = 0; col < vec.size(); col++) {
			result[row] += matrix[row][col] * vec[col];
		}
	}
	return result;
}

template<size_t precision, size_t maxIter>
double LOS<precision, maxIter>::scal(const std::vector<double>& v1, const std::vector<double>& v2) const {
	double res = 0.0;
	for (size_t ind = 0; ind < v1.size(); ind++)
		res += v1[ind] * v2[ind];
	return res;
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::zeroIter(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	
	mult(matrix, startVect, result);

	for (size_t i = 0; i < result.size(); i++) {
		r[i] = vec[i] - result[i];
		z[i] = r[i];
	}

	result = startVect;

	mult(matrix, z, p);

	curDelt = scal(r, r);
	curIter = 0;
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::iterate(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
	a = scal(p, r) / scal(p, p);
	
	for (size_t i = 0; i < result.size(); i++) {
		result[i] = result[i] + a * z[i];
		r[i] = r[i] - a * p[i];
	}
	
	mult(matrix, r, ar);
	b = -scal(p, ar) / scal(p, p);

	
	for (size_t i = 0; i < result.size(); i++) {
		z[i] = r[i] + b * z[i];
		p[i] = ar[i] + b * p[i];
	}

	mult(matrix, result, ar);
	for (size_t i = 0; i < result.size(); i++) {
		ar[i] = ar[i] - vec[i];
	}
	curDelt = scal(ar,ar);


	curIter++;
}

template<size_t precision, size_t maxIter>
bool LOS<precision, maxIter>::isEnd() {
	return (curDelt < pow(10.0, -static_cast<double>(precision))) || (curIter >= maxIter);
}

template<size_t precision, size_t maxIter>
bool LOS<precision, maxIter>::isCorrect() const {
	return curDelt < pow(10.0, -static_cast<double>(precision));
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::initStart(size_t dimention) {
	if (startVect.empty()) 
		for (size_t i = 0; i < dimention; i++)
			startVect.push_back(0.0);
	curDelt = 1E+13; // Any big number
}

template<size_t precision, size_t maxIter>
void LOS<precision, maxIter>::setStart(const std::vector<double>& start) {
	startVect = start;
}

