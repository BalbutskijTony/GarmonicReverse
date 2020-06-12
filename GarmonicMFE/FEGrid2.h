#pragma once
#include "Grid2.h"
#include "FElement2.h"
#include "MaterialMap.h"

class FEGrid2 {
public:
	FEGrid2() = default;
	~FEGrid2() = default;

	size_t size() const;
	size_t addElem(const FElement2& curElem);
	void constructFEGrid(const Grid2& geomGrid, MaterialMap* materialGrid);
	void restructMaterials(const std::vector<double>& sigmas);
	//void setParams(size_t index, const std::vector<double>& params);
	const std::vector<FElement2>& data() const;

	friend std::ostream& operator<<(std::ostream& out, const FEGrid2& grid);
private:
	std::vector<FElement2> elements;
	MaterialMap* materials;
};


