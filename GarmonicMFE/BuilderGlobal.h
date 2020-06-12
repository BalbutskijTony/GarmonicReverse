#pragma once
#include <fstream>
#include <iomanip>
#include <functional>

#include "FEGrid2.h"
#include "BasicFuncFabric.h"

#define DIM_LOCAL 4

//--------------------------------------------------------------------------------
// ядро локальной матрицы ∆Єсткости(с множителем l/6*hy/hx)
//--------------------------------------------------------------------------------
static double kernGLocalYX[DIM_LOCAL][DIM_LOCAL] = { { 2, -2, 1, -1},
											  { -2, 2, -1, 1},
											  { 1, -1, 2, -2},
											  { -1, 1, -2, 2} };

//--------------------------------------------------------------------------------
// ядро локальной матрицы ∆Єсткости(с множителем l/6*hx/hy)
//--------------------------------------------------------------------------------
static double kernGLocalXY[DIM_LOCAL][DIM_LOCAL] = { { 2, 1, -2, -1 },
											  { 1, 2, -1, -2 },
											  { -2, -1, 2, 1 },
											  { -1, -2, 1, 2 } };

//--------------------------------------------------------------------------------
// ядро локальной матрицы ћасс(ядро локальной матрицы C)
//--------------------------------------------------------------------------------
static double kernMLocal[DIM_LOCAL][DIM_LOCAL] = { { 4, 2, 2, 1 },
											{ 2, 4, 1, 2 },
											{ 2, 1, 4, 2 },
											{ 1, 2, 2, 4 } };


template<std::size_t N>
class BuilderGlobal {
public:

	void setFEGrid(const FEGrid2& newGrid);
	void setFullGrid(const Grid2& fullGrid);

	void setKernM(const double(&localKernM)[N][N]);
	void setKernG(const double(&localKernGXY)[N][N], const double(&localKernGYX)[N][N]);

	void setFSin(std::function<double(double, double, double, double, double)> fSin);
	void setFCos(std::function<double(double, double, double, double, double)> fCos);

	void initGlobalMatrix();
	void initGlobalPAndC();

	void fillGlobalMatrix();
	void fillGlobalPAndC(double w, double lambda, double hi);
	void rebuildM();
	void rebuildC(double w);

	void deinitGlobalMatrix();
	void deinitGlobalPAndC();

	void setNewParameters(const std::vector<double>& sigmas);

	std::vector<std::vector<double>> buildGlobalA() const;
	std::vector<double> buildGlobalB();

private:
	size_t getGlobalIndex(size_t feIndex, size_t localIndex) const;

	FEGrid2 feGrid;
	Grid2 fullGrid;
	double localKernM[N][N];
	double localKernGXY[N][N];
	double localKernGYX[N][N];

	double** globalM = nullptr;
	double** globalG = nullptr;
	double** globalP = nullptr;
	double** globalC = nullptr;

	double lambda, hi;

	std::function<double(double, double, double, double, double)> fSin;
	std::function<double(double, double, double, double, double)> fCos;

	BasicFuncFabric funcFabric;
	Gauss2 gauss2;
};

template<std::size_t N>
size_t BuilderGlobal<N>::getGlobalIndex(size_t feIndex, size_t localIndex) const {
	size_t countFEInRow = fullGrid.sizeX() - 1;

	size_t curFERow = feIndex / countFEInRow;
	size_t curFEColumn = feIndex % countFEInRow;

	size_t shiftByX = curFEColumn;
	size_t shiftByY = curFERow;

	size_t startGlobIndex = shiftByX + shiftByY * fullGrid.sizeX();

	return startGlobIndex + (localIndex / 2) * fullGrid.sizeX() + (localIndex % 2) * 1;
}

template<std::size_t N>
void BuilderGlobal<N>::setFEGrid(const FEGrid2& newGrid) {
	feGrid = newGrid;
}

template<std::size_t N>
void BuilderGlobal<N>::setFullGrid(const Grid2& fullGrid){
	this->fullGrid = fullGrid;
}

template<std::size_t N>
void BuilderGlobal<N>::setKernM(const double(&localKernM)[N][N]) {
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < N; col++) {
			this->localKernM[row][col] = localKernM[row][col];
		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::setKernG(const double(&localKernGXY)[N][N], const double(&localKernGYX)[N][N]) {
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < N; col++) {
			this->localKernGXY[row][col] = localKernGXY[row][col];
			this->localKernGYX[row][col] = localKernGYX[row][col];
		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::initGlobalMatrix() {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	globalM = new double* [nodeCount];
	for (size_t i = 0; i < nodeCount; i++)
		globalM[i] = new double[nodeCount]();

	globalG = new double* [nodeCount];
	for (size_t i = 0; i < nodeCount; i++)
		globalG[i] = new double[nodeCount]();
}

template<std::size_t N>
void BuilderGlobal<N>::initGlobalPAndC() {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	globalP = new double* [nodeCount];
	for (size_t i = 0; i < nodeCount; i++)
		globalP[i] = new double[nodeCount]();

	globalC = new double* [nodeCount];
	for (size_t i = 0; i < nodeCount; i++)
		globalC[i] = new double[nodeCount]();
}

template<std::size_t N>
void BuilderGlobal<N>::fillGlobalMatrix() {
	size_t feCount = feGrid.size();

	for (size_t curFE = 0; curFE < feCount; curFE++) {

		auto curElement = feGrid.data()[curFE];
		double hx = curElement.getPoint(1).x - curElement.getPoint(0).x;
		double hy = curElement.getPoint(2).y - curElement.getPoint(0).y;

		for (size_t localRow = 0; localRow < N; localRow++) {
			size_t globalRow = getGlobalIndex(curFE, localRow);

			for (size_t localCol = 0; localCol < N; localCol++) {
				size_t globalCol = getGlobalIndex(curFE, localCol);

				globalM[globalRow][globalCol] += curElement.getSigma() *
					hx * hy / 36 * kernMLocal[localRow][localCol];
				globalG[globalRow][globalCol] += /*curElement.getLambda() */ 1. / 6 * (
					hy / hx * kernGLocalYX[localRow][localCol] + hx / hy * kernGLocalXY[localRow][localCol]);
			}

		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::rebuildM() {
	size_t feCount = feGrid.size();
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	for (size_t i = 0; i < nodeCount; i++)
		for (size_t j = 0; j < nodeCount; j++)
			globalM[i][j] = 0.0;

	for (size_t curFE = 0; curFE < feCount; curFE++) {

		auto curElement = feGrid.data()[curFE];
		double hx = curElement.getPoint(1).x - curElement.getPoint(0).x;
		double hy = curElement.getPoint(2).y - curElement.getPoint(0).y;

		for (size_t localRow = 0; localRow < N; localRow++) {
			size_t globalRow = getGlobalIndex(curFE, localRow);

			for (size_t localCol = 0; localCol < N; localCol++) {
				size_t globalCol = getGlobalIndex(curFE, localCol);

				globalM[globalRow][globalCol] += curElement.getSigma() *
					hx * hy / 36 * kernMLocal[localRow][localCol];
			}

		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::rebuildC(double w) {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	for (size_t curRow = 0; curRow < nodeCount; curRow++) {
		for (size_t curCol = 0; curCol < nodeCount; curCol++) {
			globalC[curRow][curCol] = w * globalM[curRow][curCol];
		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::fillGlobalPAndC(double w, double lambda, double hi) {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	this->lambda = lambda;
	this->hi = hi;

	for (size_t curRow = 0; curRow < nodeCount; curRow++) {
		for (size_t curCol = 0; curCol < nodeCount; curCol++) {
			globalC[curRow][curCol] = w * globalM[curRow][curCol];
			globalP[curRow][curCol] = lambda * globalG[curRow][curCol];// -w * w * hi * globalM[curRow][curCol];
		}
	}
}

template<std::size_t N>
void BuilderGlobal<N>::deinitGlobalMatrix() {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	for (size_t i = 0; i < nodeCount; i++)
		delete[] globalM[i];
	delete[] globalM;

	for (size_t i = 0; i < nodeCount; i++)
		delete[] globalG[i];
	delete[] globalG;
}

template<std::size_t N>
void BuilderGlobal<N>::deinitGlobalPAndC() {
	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;

	for (size_t i = 0; i < nodeCount; i++)
		delete[] globalP[i];
	delete[] globalP;

	for (size_t i = 0; i < nodeCount; i++)
		delete[] globalC[i];
	delete[] globalC;
}

template<std::size_t N>
void BuilderGlobal<N>::setFSin(std::function<double(double, double, double, double, double)> fSin) {
	this->fSin = fSin;
}

template<std::size_t N>
void BuilderGlobal<N>::setFCos(std::function<double(double, double, double, double, double)> fCos) {
	this->fCos = fCos;
}

template<std::size_t N>
std::vector<std::vector<double>> BuilderGlobal<N>::buildGlobalA() const {
	std::vector<std::vector<double>> result;

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;
	size_t feCount = feGrid.size();

	// Resize and fill result
	{
		std::vector<double> row;
		for (size_t curNode = 0; curNode < nodeCount * 2; curNode++) {
			row.push_back(0.0);
		}
		for (size_t curRow = 0; curRow < nodeCount * 2; curRow++) {
			result.push_back(row);
		}
	}

	for (size_t curRow = 0; curRow < nodeCount; curRow++) {
		for (size_t curCol = 0; curCol < nodeCount; curCol++) {
			result[2 * curRow][2 * curCol] += globalP[curRow][curCol];
			result[2 * curRow + 1][2 * curCol + 1] += globalP[curRow][curCol];
			result[2 * curRow][2 * curCol + 1] += -globalC[curRow][curCol];
			result[2 * curRow + 1][2 * curCol] += globalC[curRow][curCol];
		}
	}

	std::ofstream out("globalC.txt");
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			out << globalC[i][j] << '\t';
		}
		out << '\n';
	}
	out.close();

	out.open("globalP.txt");
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			out << globalP[i][j] << '\t';
		}
		out << '\n';
	}
	out.close();

	out.open("globalG.txt");
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			out << globalG[i][j] << '\t';
		}
		out << '\n';
	}
	out.close();

	out.open("globalM.txt");
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			out << globalM[i][j] << '\t';
		}
		out << '\n';
	}
	out.close();

	return result;
}

template<std::size_t N>
std::vector<double> BuilderGlobal<N>::buildGlobalB() {
	std::vector<double> result;

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();
	size_t nodeCount = countNodeX * countNodeY;
	size_t feCount = feGrid.size();

	// Resize and fill result
	for (size_t curNode = 0; curNode < nodeCount * 2; curNode++) {
		result.push_back(0.0);
	}

	for (size_t curFE = 0; curFE < feCount; curFE++) {

		auto curFEData = feGrid.data()[curFE];
		double hx = curFEData.getPoint(1).x - curFEData.getPoint(0).x;
		double hy = curFEData.getPoint(2).y - curFEData.getPoint(0).y;

		gauss2.setArea(curFEData.getPoint(0).x, curFEData.getPoint(1).x, curFEData.getPoint(0).y, curFEData.getPoint(2).y);
		funcFabric.setArea(curFEData.getPoint(0).x, curFEData.getPoint(1).x, curFEData.getPoint(0).y, curFEData.getPoint(2).y);

		// TODO: ѕрверить индексы
		for (size_t local = 0; local < N; local++) {
			size_t startGlobal = getGlobalIndex(curFE, local);
			std::function<double(double, double)> fSinLoc = [=](double x, double y)
			{return fSin(x, y, feGrid.data()[curFE].getSigma(), lambda, hi); };
			std::function<double(double, double)> fCosLoc = [=](double x, double y)
			{return fCos(x, y, feGrid.data()[curFE].getSigma(), lambda, hi); };

			result[2 * startGlobal] += gauss2.integrate(fSinLoc, funcFabric.getBasic(local));
				
			result[2 * startGlobal + 1] += gauss2.integrate(fCosLoc, funcFabric.getBasic(local));

			// TODO: ќшибка скорее всего тут
			//for (size_t localCol = 0; localCol < N; localCol++) {
				//result[2 * startGlobal] += fSin(curFEData.getPoint(localCol).x, curFEData.getPoint(localCol).y)
				//	* globalC[startGlobal][getGlobalIndex(curFE, localCol)];
				//result[2 * startGlobal + 1] += fCos(curFEData.getPoint(localCol).x, curFEData.getPoint(localCol).y)
				//	* globalC[startGlobal][getGlobalIndex(curFE, localCol)];

			//}
		}
	}

	return result;
}

template<std::size_t N>
void BuilderGlobal<N>::setNewParameters(const std::vector<double>& sigmas) {

	feGrid.restructMaterials(sigmas);
}