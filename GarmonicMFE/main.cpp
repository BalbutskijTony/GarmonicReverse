
#include <iostream>
//
//#include "FEGrid2.h"
//
//#include "LOS.h"
//#include "StabilizedBiorderedGradient.h"
//
//#include "Gauss2.h"
//#include "BuilderGlobal.h"

#include "GarmonicTask.h"
#include "ReverseSolver.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "MonolithMap.h"
#include "LayeredMap.h"

// DEPRECATED START
// Метод Гаусса
std::vector<double> gauss(const std::vector<std::vector<double>>& matrix, const std::vector<double>& right) {
	auto a = matrix;
	auto b = right;
	std::vector<double> result(b.size());

	int i, k, m, im; long double amm, aim, r;

	for (m = 0; m < b.size() - 1; m++) {// m
		amm = a[m][m]; im = m;

		for (i = m; i < b.size(); i++)
			if (a[i][m] > amm) { amm = a[i][m]; im = i; }

		if (im != m) {
			r = b[im]; b[im] = b[m]; b[m] = r;
			for (k = m; k < b.size(); k++)
			{
				r = a[im][k]; a[im][k] = a[m][k]; a[m][k] = r;
			}
		}

		for (k = m; k < b.size(); k++)
			a[m][k] = a[m][k] / amm; // 3.16

		b[m] = b[m] / amm; //

		for (i = m + 1; i < b.size(); i++) {// i
			aim = a[i][m];
			for (k = m; k < b.size(); k++)
				a[i][k] = a[i][k] - a[m][k] * aim; // 3.17
			b[i] = b[i] - b[m] * aim; //
		}// end i
	}// end m

	result[b.size() - 1] = b[b.size() - 1] / a[b.size() - 1][b.size() - 1]; // 3.19
	for (i = b.size() - 1; i >= 0; i--) {// i
		result[i] = b[i]; //
		for (k = i + 1; k < b.size(); k++) // 3.20
			result[i] = result[i] - a[i][k] * result[k]; //
	}// end i

	return result;
}

std::vector<double> solveSLAE(const std::vector<std::vector<double>> matrix, const std::vector<double> rightPart) {
	auto result = rightPart;
	auto _matrix = matrix;
	double tmp;

	for (int i = 0; i < result.size(); i++) {
		tmp = _matrix[i][i];
		for (int j = i; j < result.size(); j++)
			_matrix[i][j] /= tmp;
		result[i] /= tmp;
		for (int j = i + 1; j < result.size(); j++) {
			tmp = _matrix[j][i];
			for (int k = i; k < result.size(); k++)
				_matrix[j][k] -= tmp * _matrix[i][k];
			result[j] -= tmp * result[i];
		}
		for (int j = i - 1; j >= 0; j--) {
			tmp = _matrix[j][i];
			for (int k = i; k < result.size(); k++)
				_matrix[j][k] -= tmp * _matrix[i][k];
			result[j] -= tmp * result[i];
		}
	}

	return result;
}


void applyLeftCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart,
	std::function<double(double, double)> valueSin, std::function<double(double, double)> valueCos) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curLayer = 0; curLayer < countNodeY; curLayer++) {
		
		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curLayer * countNodeX][curCol] = 0.0;
			matrix[2 * curLayer * countNodeX + 1][curCol] = 0.0;
		}

		matrix[2 * curLayer * countNodeX][2 * curLayer * countNodeX] = 1.0;// *1.E7;
		matrix[2 * curLayer * countNodeX + 1][2 * curLayer * countNodeX + 1] = 1.0;// *1.E7;

		rightPart[2 * curLayer * countNodeX] = valueSin(fullGrid.getNode(0, curLayer).first, fullGrid.getNode(0, curLayer).second);// *1.E7;
		rightPart[2 * curLayer * countNodeX + 1] = valueCos(fullGrid.getNode(0, curLayer).first, fullGrid.getNode(0, curLayer).second);// *1.E7;

	}

}

void applyRightCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart,
	std::function<double(double, double)> valueSin, std::function<double(double, double)> valueCos) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curLayer = 0; curLayer < countNodeY; curLayer++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX][curCol] = 0.0;
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1][curCol] = 0.0;
		}

		matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX]
			[2 * (countNodeX - 1) + 2 * curLayer * countNodeX] = 1.0;// *1.E7;
		matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1]
			[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1] = 1.0;// *1.E7;

		rightPart[2 * (countNodeX - 1) + 2 * curLayer * countNodeX] =
			valueSin(fullGrid.getNode(countNodeX - 1, curLayer).first, fullGrid.getNode(countNodeX - 1, curLayer).second);// *1.E7;
		rightPart[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1] =
			valueCos(fullGrid.getNode(countNodeX - 1, curLayer).first, fullGrid.getNode(countNodeX - 1, curLayer).second);// *1.E7;

	}
}

void applyBottomCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart,
	std::function<double(double, double)> valueSin, std::function<double(double, double)> valueCos) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curX = 0; curX < countNodeX; curX++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curX][curCol] = 0.0;
			matrix[2 * curX + 1][curCol] = 0.0;
		}

		matrix[2 * curX][2 * curX] = 1.0;// *1.E7;
		matrix[2 * curX + 1][2 * curX + 1] = 1.0;// *1.E7;

		rightPart[2 * curX] = valueSin(fullGrid.getNode(curX, 0).first, fullGrid.getNode(curX, 0).second);// *1.E7;
		rightPart[2 * curX + 1] = valueCos(fullGrid.getNode(curX, 0).first, fullGrid.getNode(curX, 0).second);// *1.E7;

	}
}

void applyTopCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart,
	std::function<double(double, double)> valueSin, std::function<double(double, double)> valueCos) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curX = 0; curX < countNodeX; curX++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][curCol] = 0.0;
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][curCol] = 0.0;
		}

		matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][2 * curX + 2 * countNodeX * (countNodeY - 1)] = 1.0;// *1.E7;
		matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][2 * curX + 2 * countNodeX * (countNodeY - 1) + 1] = 1.0;// *1.E7;

		rightPart[2 * curX + 2 * countNodeX * (countNodeY - 1)] =
			valueSin(fullGrid.getNode(curX, countNodeY - 1).first, fullGrid.getNode(curX, countNodeY - 1).second);// *1.E7;
		rightPart[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1] =
			valueCos(fullGrid.getNode(curX, countNodeY - 1).first, fullGrid.getNode(curX, countNodeY - 1).second);// *1.E7;

	}
}

// DEPRECATED END

// U sin real = xy
// U cos real = xy
// f sin real = - w sigma x y - w^2 hi x y
// f cos real = w sigma x y - w^2 hi x y 
//void generateSyntetic1(BuilderGlobal<4>& builder, double w, double lambda, double hi, double sigma) {
//
//	builder.setFCos([=](double x, double y) {return w * sigma * x * y - w * w * hi * x * y; });
//	builder.setFSin([=](double x, double y) {return - w * sigma * x * y - w * w * hi * x * y; });
//}


int main() {

	// DEPRECATED
	Grid2 mainGrid;
	mainGrid.fromFile("mainGrid.txt");
	int countByX = mainGrid.sizeX() - 1;
	int countByY = mainGrid.sizeY() - 1;

	Grid2 globalGrid = mainGrid.procreate("xKoeff.txt", "yKoeff.txt");

	//FEGrid2 feGrid;
	//feGrid.constructFEGrid(globalGrid);

	//double w = 1;
	//double lambda = 0;
	//double hi = 0;
	//double sigma = 1;

	//BuilderGlobal<DIM_LOCAL> builder;

	//builder.setFEGrid(feGrid);
	//builder.setFullGrid(globalGrid);
	//builder.setKernM(kernMLocal);
	//builder.setKernG(kernGLocalXY, kernGLocalYX);

	//builder.initGlobalMatrix();
	//builder.initGlobalPAndC();

	//builder.fillGlobalMatrix();
	//builder.fillGlobalPAndC(w, lambda, hi, sigma);

	//generateSyntetic1(builder, w, lambda, hi, sigma);
	double w = 2. * M_PI * 0.01;
	double hi = 0.0;
	double lambda = 1.E+7 / (4. * M_PI);
	
	//MaterialMap* matMap = new MonolithMap();
	MaterialMap* matMap = new LayeredMap(); matMap->setColCount(1);
	matMap->loadFromFile("material.txt");

	std::function<double(double, double, double, double, double)> fSin = [=](double x, double y, double sigma, double lambda, double hi)
	{return 0.0; };
	std::function<double(double, double, double, double, double)> fCos = [=](double x, double y, double sigma, double lambda, double hi)
	{return 0.0; };

	std::function<double(double, double)> USin =
		[](double x, double y) { return x * y; };
	std::function<double(double, double)> UCos =
		[](double x, double y) { return x * y; };

	// TEST 2
	//std::function<double(double, double, double, double, double)> fSin = [=](double x, double y, double sigma, double lambda, double hi)
	//{return -w * sigma * (-3. * x * y - 4 * x) - w * w * hi * (x * y + 2 * y); };
	//std::function<double(double, double, double, double, double)> fCos = [=](double x, double y, double sigma, double lambda, double hi)
	//{return w * sigma * (x * y + 2 * y) - w * w * hi * (-3. * x * y - 4 * x); };

	//std::function<double(double, double)> USin =
	//	[](double x, double y) { return x * y + 2 * y; };
	//std::function<double(double, double)> UCos =
	//	[](double x, double y) { return -3. * x * y - 4 * x; };

	// TEST 3
	//std::function<double(double, double, double, double, double)> fSin = [=](double x, double y, double sigma, double lambda, double hi)
	//{return -6. * lambda * (x * y * y * y + y * x * x * x)//-2. * lambda * (x * x + y * y)//-2. * lambda * (y * y * y + 1. + 3. * y * x * x)//
	//	- w * sigma * x * x * x * y * y * y//- w * sigma * (-x * x * y * y + 4.)//-w * sigma * 2 * (x * x * x + 4) * (y * y - 3) // 
	//	- w * w * hi * x * x * x * y * y * y; };//- w * w * hi * x * x * y * y; };//- w * w * hi * x * x * (y * y * y + 1); }; //
	//std::function<double(double, double, double, double, double)> fCos = [=](double x, double y, double sigma, double lambda, double hi)
	//{return -6. * lambda * (x * y * y * y + y * x * x * x)//2. * lambda * (x * x + y * y)//-2. * lambda * (6. * x * (y*y - 3) + 2 *(x*x*x + 4)) //
	//	+ w * sigma * x * x * x * y * y * y //+ w * sigma * (x * x * y * y)//+ w * sigma * x * x * (y * y * y + 1) //
	//	- w * w * hi * x * x * x * y * y * y; }; //- w * w * hi * (4. - x * x * y * y); };//- w * w * hi * 2 * (x * x * x + 4) * (y * y - 3); }; // 

	//std::function<double(double, double)> USin =
	//	[](double x, double y) { return x * x * x * y * y * y; };//return x * x * y * y; };//x * x * (y * y * y  + 1); }; //
	//std::function<double(double, double)> UCos =
	//	[](double x, double y) { return x * x * x * y * y * y; };//return -x * x * y * y + 4.; };//2 * (x*x*x + 4) * (y*y - 3); }; //

 	GarmonicTask testTask;
	testTask.setMaterialType(GarmonicTask::MaterialType::globalSigma);
	testTask.createGrid("mainGrid.txt", "xKoeff.txt", "yKoeff.txt", matMap); //???
	testTask.setRightFunction(fSin, fCos);
	testTask.setHi(hi);
	testTask.setLambda(lambda);
	testTask.setW(w);
	testTask.init();
	//testTask.setParameters(matMap.data()); //???
	testTask.setBoundaryCondition(
		[](double x, double y) {return 1; }, GarmonicTask::ConditionType::condition2nd,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::condition2nd,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::condition1st,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::condition1st,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::noCondition, // ::condition1st,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::condition1st,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::noCondition,
		[](double x, double y) {return 0; }, GarmonicTask::ConditionType::condition1st);
	testTask.solve();
	auto solution = testTask.getSolution();
	auto realRecievers = testTask.getReceiveValues();

	//std::cout << "Correct " << std::boolalpha << testTask.isCorrect() <<std::endl;

	ReverseSolver revSolver;
	revSolver.setStartParams({ 
		0.001,
		0.001,
		0.001,
		0.001/*,
		0.0009,
		0.0009*/
		/*0.000001,
		0.000001,
		0.000001*/
		//0.01,
		//0.01,
		//0.01,
		//0.01,
		//0.01,//
		//0.01,//
		//0.01,//
		//0.01,//
		//0.01,
		//0.01,
		//0.01,
		//0.01 
	}); // , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
	revSolver.setError(1.E-7);
	revSolver.setMaxIter(10000);
	revSolver.setTask(&testTask);
	revSolver.setTrueEpsilon(testTask.getReceiveValues());
	revSolver.startState();
	revSolver.solve();
	auto param = revSolver.getActualParameters();


	std::ofstream res("result.txt");

	res << "X\t\tY\t\tUsc\t\tUcc\t\tUsr\t\tUcr" << std::endl;
	res << std::fixed << std::setprecision(10);
	for (int y = 0; y < globalGrid.sizeY(); y++) {
		for (int x = 0; x < globalGrid.sizeX(); x++) {
			res << globalGrid.getNode(x, y).first << '\t' <<
				globalGrid.getNode(x, y).second << '\t' <<
				solution[2 * x + 2 * globalGrid.sizeX() * y] << '\t' <<
				solution[2 * x + 2 * globalGrid.sizeX() * y + 1] << std::endl;
		}
	}

	res.close();

	system("pause");
	return 0;
}
