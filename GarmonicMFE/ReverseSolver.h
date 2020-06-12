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
	// Выставление истинного значения функционала
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

	// Текущее значение решения
	std::vector<double> curSolve;
	// Истинное решение
	std::vector<double> originSolve;
	// Вектор для хранения производных Решения
	std::vector<std::vector<double>> shiftSolve;
	// Хранит дельта-P по всем координатам
	std::vector<double> shiftParams;
	// Текущий вектор параметров
	std::vector<double> curParams;
	// Вектор весов для влияния компонент решения
	std::vector<double> weight;

	double calcError(const std::vector<double>& solution) const;

	void outData(const std::vector<std::vector<double>>& matrix, const std::string& fileName) const;
	void outData(const std::vector<double>& vector, const std::string& fileName) const;
};

