#include "GarmonicTask.h"

// Task's methods start
void GarmonicTask::init() {
	builder.setFEGrid(feGrid);
	builder.setFullGrid(fullGrid);
	builder.setKernM(kernMLocal);
	builder.setKernG(kernGLocalXY, kernGLocalYX);

	builder.initGlobalMatrix();
	builder.initGlobalPAndC();

	builder.fillGlobalMatrix();
}

void GarmonicTask::setParameters(const std::vector<double>& params) {
	// TODO: Временное решение

	builder.setNewParameters(params);
	builder.rebuildM();
	builder.rebuildC(w);
}

void GarmonicTask::solve() {

	builder.fillGlobalPAndC(w, lambda, hi);
	globalA = builder.buildGlobalA();
	globalB = builder.buildGlobalB();

	std::ofstream out("MatrixA.txt");
	for (size_t i = 0; i < globalA.size(); i++)
	{
		for (size_t j = 0; j < globalA.size(); j++) {
			out << globalA[i][j] << '\t';
		}
		out << std::endl;
	}
	out.close();

	std::ofstream outVect("VectorB.txt");
	for (size_t i = 0; i < globalA.size(); i++)
	{
		out << globalB[i] << '\t';
	}
	outVect.close();

	if(ConditionType::noCondition != topSinType)
		applyTopCondition(fullGrid, globalA, globalB);
	if (ConditionType::noCondition != bottomSinType)
		applyBottomCondition(fullGrid, globalA, globalB);
	if (ConditionType::noCondition != leftSinType)
		applyLeftCondition(fullGrid, globalA, globalB);
	if (ConditionType::noCondition != rightSinType)
		applyRightCondition(fullGrid, globalA, globalB);

	//solver.setEndType(StabilizedBiorderedGradient::TypeEnd::BY_ERROR);

	Eigen::MatrixXd matrix(globalA.size(), globalA.size());
	Eigen::VectorXd vect(globalA.size());

	for (size_t row = 0; row < globalA.size(); row++) {
		for (size_t col = 0; col < globalA.size(); col++) {
			matrix(row, col) = globalA[row][col];
		}

		vect(row) = globalB[row];
	}

	Eigen::MatrixXd tmp = matrix.lu().solve(vect);
	Eigen::MatrixXd tmp1 = matrix.householderQr().solve(vect);
	result.clear();
	for (size_t row = 0; row < globalB.size(); row++)
		result.push_back(tmp(row, 0));
	//result = solver.solve(globalA, globalB);
}

const std::vector<double>& GarmonicTask::getSolution() const {
	return result;
}

const std::vector<double>& GarmonicTask::getReceiveValues() {
	receive.clear();
	// Вытаскиваем верхний слой узлов сетки без крайних элементов
	// и умножаем на -omega
	for (size_t elInd = result.size() - 2 * fullGrid.sizeX(); elInd < result.size(); elInd++) {
	//for (size_t elInd = result.size() - 2 * fullGrid.sizeX() + 2; elInd < result.size() - 2 * fullGrid.sizeX() + 3; elInd++) {
		receive.push_back(result[elInd]); // * -w
	}
	return receive;
}

size_t GarmonicTask::paramCount() const {
	if (MaterialType::globalSigma == materialType) return 1;
	else return feGrid.size();
}
// Task's methods end

void GarmonicTask::setRightFunction(const std::function<double(double, double, double, double, double)>& fSin, 
	const std::function<double(double, double, double, double, double)>& fCos) {
	builder.setFSin(fSin);
	builder.setFCos(fCos);
}

void GarmonicTask::createGrid(const std::string& gridFile, const std::string& xKoeff, const std::string& yKoeff, MaterialMap* map) {
	mainGrid.fromFile(gridFile);
	fullGrid = mainGrid.procreate(xKoeff, yKoeff);

	feGrid.constructFEGrid(fullGrid, map);


	builder.setFEGrid(feGrid);
	builder.setFullGrid(fullGrid);

	builder.initGlobalMatrix();
	builder.initGlobalPAndC();
}

void GarmonicTask::changeMaterialGrid(MaterialMap* map) {

	feGrid.constructFEGrid(fullGrid, map);
	
}

void GarmonicTask::setBoundaryCondition(const std::function<double(double, double)>& topSin, ConditionType topSinType,
	const std::function<double(double, double)>& topCos, ConditionType topCosType,
	const std::function<double(double, double)>& bottomSin, ConditionType bottomSinType,
	const std::function<double(double, double)>& bottomCos, ConditionType bottomCosType,
	const std::function<double(double, double)>& leftSin, ConditionType leftSinType,
	const std::function<double(double, double)>& leftCos, ConditionType leftCosType,
	const std::function<double(double, double)>& rightSin, ConditionType rightSinType,
	const std::function<double(double, double)>& rightCos, ConditionType rightCosType) {

	this->topConditionSin = topSin;
	this->topConditionCos = topCos;
	this->bottomConditionSin = bottomSin;
	this->bottomConditionCos = bottomCos;
	this->leftConditionSin = leftSin;
	this->leftConditionCos = leftCos;
	this->rightConditionSin = rightSin;
	this->rightConditionCos = rightCos;

	this->topSinType = topSinType;
	this->topCosType = topCosType;
	this->bottomSinType = bottomSinType;
	this->bottomCosType = bottomCosType;
	this->leftSinType = leftSinType;
	this->leftCosType = leftCosType;
	this->rightSinType = rightSinType;
	this->rightCosType = rightCosType;
}

// Temporary functions
void GarmonicTask::setLambda(double lambda) {
	this->lambda = lambda;
}

void GarmonicTask::setHi(double hi) {
	this->hi = hi;
}

void GarmonicTask::setW(double w) {
	this->w = w;
}


void GarmonicTask::applyLeftCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curLayer = 0; curLayer < countNodeY; curLayer++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curLayer * countNodeX][curCol] = 0.0;
			matrix[2 * curLayer * countNodeX + 1][curCol] = 0.0;
		}

		if (ConditionType::condition1st == leftSinType) {
			matrix[2 * curLayer * countNodeX][2 * curLayer * countNodeX] = 1.0;
		}
		else {
			double stepX = fullGrid.getNode(1, 0).first - fullGrid.getNode(0, 0).first; // TODO: Слабое место для расширения
			matrix[2 * curLayer * countNodeX][2 * curLayer * countNodeX] = -lambda / stepX;
			matrix[2 * curLayer * countNodeX][2 * curLayer * countNodeX + 2] = lambda / stepX;
		}
		rightPart[2 * curLayer * countNodeX] = leftConditionSin(fullGrid.getNode(0, curLayer).first, fullGrid.getNode(0, curLayer).second);

		if (ConditionType::condition1st == leftCosType) {
			matrix[2 * curLayer * countNodeX + 1][2 * curLayer * countNodeX + 1] = 1.0;// *1.E7;
		}
		else {
			double stepX = fullGrid.getNode(1, 0).first - fullGrid.getNode(0, 0).first; // TODO: Слабое место для расширения
			matrix[2 * curLayer * countNodeX + 1][2 * curLayer * countNodeX + 1] = -lambda / stepX;
			matrix[2 * curLayer * countNodeX + 1][2 * curLayer * countNodeX + 3] = lambda / stepX;
		}
		rightPart[2 * curLayer * countNodeX + 1] = leftConditionCos(fullGrid.getNode(0, curLayer).first, fullGrid.getNode(0, curLayer).second);// *1.E7;

	}

}

void GarmonicTask::applyRightCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curLayer = 0; curLayer < countNodeY; curLayer++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX][curCol] = 0.0;
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1][curCol] = 0.0;
		}

		if (ConditionType::condition1st == rightSinType) {
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX] = 1.0;
		}
		else {
			double stepX = fullGrid.getNode(countNodeX - 1, 0).first - fullGrid.getNode(countNodeX - 2, 0).first; // TODO: Слабое место для расширения
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX] = lambda / stepX;
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX - 2] = -lambda / stepX;
		}
		rightPart[2 * (countNodeX - 1) + 2 * curLayer * countNodeX] =
			rightConditionSin(fullGrid.getNode(countNodeX - 1, curLayer).first, fullGrid.getNode(countNodeX - 1, curLayer).second);

		if (ConditionType::condition1st == rightSinType) {
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1] = 1.0;// *1.E7;
		}
		else {
			double stepX = fullGrid.getNode(countNodeX - 1, 0).first - fullGrid.getNode(countNodeX - 2, 0).first; // TODO: Слабое место для расширения
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1] = lambda / stepX;
			matrix[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1]
				[2 * (countNodeX - 1) + 2 * curLayer * countNodeX - 1] = -lambda / stepX;
		}
		rightPart[2 * (countNodeX - 1) + 2 * curLayer * countNodeX + 1] =
			rightConditionCos(fullGrid.getNode(countNodeX - 1, curLayer).first, fullGrid.getNode(countNodeX - 1, curLayer).second);// *1.E7;

	}
}

void GarmonicTask::applyBottomCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curX = 0; curX < countNodeX; curX++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curX][curCol] = 0.0;
			matrix[2 * curX + 1][curCol] = 0.0;
		}

		if (ConditionType::condition1st == bottomSinType) {
			matrix[2 * curX][2 * curX] = 1.0;
		}
		else {
			double stepY = fullGrid.getNode(0, 1).second - fullGrid.getNode(0, 0).second; // TODO: Слабое место для расширения
			matrix[2 * curX][2 * curX] = -lambda / stepY;
			matrix[2 * curX][2 * curX + 2 * countNodeX] = lambda / stepY;
		}
		rightPart[2 * curX] = bottomConditionSin(fullGrid.getNode(curX, 0).first, fullGrid.getNode(curX, 0).second);

		if (ConditionType::condition1st == bottomCosType) {
			matrix[2 * curX + 1][2 * curX + 1] = 1.0;// *1.E7;
		}
		else {
			double stepY = fullGrid.getNode(0, 1).second - fullGrid.getNode(0, 0).second; // TODO: Слабое место для расширения
			matrix[2 * curX + 1][2 * curX + 1] = -lambda / stepY;
			matrix[2 * curX + 1][2 * curX + 2 * countNodeX + 1] = lambda / stepY;
		}
		rightPart[2 * curX + 1] = bottomConditionCos(fullGrid.getNode(curX, 0).first, fullGrid.getNode(curX, 0).second);

	}
}

void GarmonicTask::applyTopCondition(const Grid2& fullGrid,
	std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart) {

	size_t countNodeX = fullGrid.sizeX();
	size_t countNodeY = fullGrid.sizeY();

	for (size_t curX = 0; curX < countNodeX; curX++) {

		for (size_t curCol = 0; curCol < matrix.size(); curCol++) {
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][curCol] = 0.0;
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][curCol] = 0.0;
		}

		if (ConditionType::condition1st == topCosType) {
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][2 * curX + 2 * countNodeX * (countNodeY - 1)] = 1.0;
		}
		else {
			double stepY = fullGrid.getNode(0, countNodeY - 1).second - fullGrid.getNode(0, countNodeY - 2).second; // TODO: Слабое место для расширения
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][2 * curX + 2 * countNodeX * (countNodeY - 1)] = lambda / stepY;
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1)][2 * curX + 2 * countNodeX * (countNodeY - 2)] = -lambda / stepY;
		}
		rightPart[2 * curX + 2 * countNodeX * (countNodeY - 1)] =
			topConditionSin(fullGrid.getNode(curX, countNodeY - 1).first, fullGrid.getNode(curX, countNodeY - 1).second);// *1.E7;

		if (ConditionType::condition1st == topSinType) {
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][2 * curX + 2 * countNodeX * (countNodeY - 1) + 1] = 1.0;
		}
		else {
			double stepY = fullGrid.getNode(0, countNodeY - 1).second - fullGrid.getNode(0, countNodeY - 2).second; // TODO: Слабое место для расширения
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][2 * curX + 2 * countNodeX * (countNodeY - 1) + 1] = lambda / stepY;
			matrix[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1][2 * curX + 2 * countNodeX * (countNodeY - 2) + 1] = -lambda / stepY;
		}
		rightPart[2 * curX + 2 * countNodeX * (countNodeY - 1) + 1] =
			topConditionCos(fullGrid.getNode(curX, countNodeY - 1).first, fullGrid.getNode(curX, countNodeY - 1).second);

	}
}

bool GarmonicTask::isCorrect() const {
	return true;// solver.isCorrect();
}

void GarmonicTask::setMaterialType(MaterialType newType) {
	materialType = newType;
}