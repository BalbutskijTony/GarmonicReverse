#pragma once
#include <vector>
#include <functional>

class StabilizedBiorderedGradient
{
public:
	StabilizedBiorderedGradient() = default;
	~StabilizedBiorderedGradient() = default;

	void setStart(const std::vector<double>& start);
	const std::vector<double>& solve(const std::vector<std::vector<double>>& matrix, const std::vector<double>& right);
	
	void setMaxIter(size_t iterCount);
	void setCritW(double wCrit);
	void setCritRNorm(double norm);
	
	bool isCorrect() const;

	size_t getIterCount() const;

	enum class TypeEnd {
		BY_OMEGA,
		BY_ERROR,
		ONLY_BY_ITER
	};

	void setEndType(TypeEnd type);
private:
	//----------------------------------------------------------------------------------------------------------------------
	void initAlocate(size_t dimention);
	void iterate(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
	bool isEnd();
	void zeroIter(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
	//----------------------------------------------------------------------------------------------------------------------
	void copy(std::vector<double>& to, const std::vector<double>& from);
	void mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec, std::vector<double>& result);
	std::vector<double> mult(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec);
	double scal(const std::vector<double>& v1, const std::vector<double>& v2) const;
	
	double norm(const std::vector<double>& v);
	//----------------------------------------------------------------------------------------------------------------------
	void initZero(std::vector<double>& vect, size_t count);

	size_t maxIter = 0;
	size_t curIter = 0;
	std::vector<double> r, rTild, p, v, t, s, x;
	double a, w, roCur, roLast, b;
	double criticalW;
	double criticalE;

	const std::vector<double>* right;

	TypeEnd endType;
	std::function<bool()> endCondition;

};

