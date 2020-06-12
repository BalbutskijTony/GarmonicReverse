#pragma once
#include "Task.h"

//#include "StabilizedBiorderedGradient.h"

#include "LOS.h"

#include "Gauss2.h"

#include "BuilderGlobal.h"
#include "FEGrid2.h"

#include "Eigen/Dense"

class GarmonicTask : public Task
{
public:
	enum class MaterialType {
		globalSigma,
		localSigma
	};

	enum class ConditionType {
		condition1st,
		condition2nd,
		noCondition
	};

	GarmonicTask() : materialType(MaterialType::localSigma) {};
	virtual ~GarmonicTask() override = default;

	virtual void init() override;
	virtual void setParameters(const std::vector<double>& params) override;
	virtual void solve() override;
	virtual const std::vector<double>& getReceiveValues() override;
	const std::vector<double>& getSolution() const;
	virtual size_t paramCount() const override;

	void setRightFunction(const std::function<double(double, double, double, double, double)>& fSin,
		const std::function<double(double, double, double, double, double)>& fCos);
	void createGrid(const std::string& gridFile, const std::string& xKoeff, const std::string& yKoeff, MaterialMap* map);
	void changeMaterialGrid(MaterialMap* map);

	void setBoundaryCondition(const std::function<double(double, double)>& topSin, ConditionType topSinType,
		const std::function<double(double, double)>& topCos, ConditionType topCosType,
		const std::function<double(double, double)>& bottomSin, ConditionType bottomSinType,
		const std::function<double(double, double)>& bottomCos, ConditionType bottomCosType,
		const std::function<double(double, double)>& leftSin, ConditionType leftSinType,
		const std::function<double(double, double)>& leftCos, ConditionType leftCosType,
		const std::function<double(double, double)>& rightSin, ConditionType rightSinType,
		const std::function<double(double, double)>& rightCos, ConditionType rightCosType);

	// TODO: temporary functions start
	void setLambda(double lambda);
	void setHi(double hi);
	void setW(double w);

	bool isCorrect() const;

	void setMaterialType(MaterialType newType);

private:

	// Functions region START
	std::function<double(double, double)> topConditionSin;
	std::function<double(double, double)> topConditionCos;
	std::function<double(double, double)> bottomConditionSin;
	std::function<double(double, double)> bottomConditionCos;
	std::function<double(double, double)> leftConditionSin;
	std::function<double(double, double)> leftConditionCos;
	std::function<double(double, double)> rightConditionSin;
	std::function<double(double, double)> rightConditionCos;

	ConditionType topSinType;
	ConditionType topCosType;
	ConditionType bottomSinType;
	ConditionType bottomCosType;
	ConditionType leftSinType;
	ConditionType leftCosType;
	ConditionType rightSinType;
	ConditionType rightCosType;
	// Functions region END
	Grid2 mainGrid;
	Grid2 fullGrid;
	double w, lambda, hi, sigma;
	FEGrid2 feGrid;
	BuilderGlobal<DIM_LOCAL> builder;
	//StabilizedBiorderedGradient solver;

	MaterialType materialType;

	std::vector<std::vector<double>> globalA;
	std::vector<double> globalB;
	std::vector<double> result;

	std::vector<double> receive;

	void applyLeftCondition(const Grid2& fullGrid,
		std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart);
	void applyRightCondition(const Grid2& fullGrid,
		std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart);
	void applyBottomCondition(const Grid2& fullGrid,
		std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart);
	void applyTopCondition(const Grid2& fullGrid,
		std::vector<std::vector<double>>& matrix, std::vector<double>& rightPart);
};

