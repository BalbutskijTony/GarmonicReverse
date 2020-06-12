#pragma once
#include "Task.h"
#include "StabilizedBiorderedGradient.h"
#include "LOS.h"
#include <algorithm>

#include <Eigen/Dense>
//#include "..\Eigen\src\Cholesky\LLT.h"
//#include "..\Eigen\src\Core\Matrix.h"

#include <fstream>



class ReverseSolver
{
public:
	void setTask(Task* newTask);
	void setStartParams(const std::vector<double>& params);
	void solve();
	void oldSolve();
	std::vector<double> getActualParameters();

	void startState();
	// ����������� ��������� �������� �����������
	void setTrueEpsilon(const std::vector<double>& epsilon);
	void setError(double error);

	void setMaxIter(size_t maxIter);
private:
	void iterate();
	bool isEnd();

	LOS<20,10000> solver;
	//StabilizedBiorderedGradient solver;
	

	size_t curIter;
	size_t maxIter = 100;

	double error = 1.E-7;
	double curError;

	Task* curTask;
	std::vector<std::vector<double>> jacobi;
	std::vector<double> right;

	// ������� �������� �������
	std::vector<double> curSolve;
	// �������� �������
	std::vector<double> originSolve;
	// ������ ��� �������� ����������� �������
	std::vector<std::vector<double>> shiftSolve;
	// ������ ������-P �� ���� �����������
	std::vector<double> shiftParams;
	// ������� ������ ����������
	std::vector<double> curParams;
	// ������ ����� ��� ������� ��������� �������
	std::vector<double> weight;

	double calcError(const std::vector<double>& solution) const;

	void outData(const std::vector<std::vector<double>>& matrix, const std::string& fileName) const;
	void outData(const std::vector<double>& vector, const std::string& fileName) const;
};

