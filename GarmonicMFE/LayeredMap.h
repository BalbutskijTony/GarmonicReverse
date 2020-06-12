#pragma once
#include "MaterialMap.h"

class LayeredMap : public MaterialMap
{
public:
	virtual ~LayeredMap() override = default;
	virtual void loadFromFile(const std::string& filename) override;

	virtual void setRowCount(size_t count) override;
	virtual void setColCount(size_t count) override;

	virtual double getSigma(size_t index) const override;
	virtual double getSigma(size_t row, size_t col) const override;

	virtual void setSigma(size_t index, double sigma) override;
	virtual void setSigma(size_t row, size_t col, double sigma) override;

	virtual void restructMaterial(const std::vector<double>& materials) override;
};

