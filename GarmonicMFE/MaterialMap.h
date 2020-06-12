#pragma once
#include <fstream>
#include <string>
#include <vector>

class MaterialMap abstract
{
public:
	MaterialMap() = default;
	virtual ~MaterialMap() = default;
	virtual void loadFromFile(const std::string& filename) = 0;
	
	virtual void setRowCount(size_t count) = 0;
	virtual void setColCount(size_t count) = 0;

	virtual double getSigma(size_t index) const = 0;
	virtual double getSigma(size_t row, size_t col) const = 0;

	virtual void setSigma(size_t index, double sigma) = 0;
	virtual void setSigma(size_t row, size_t col, double sigma) = 0;

	virtual void restructMaterial(const std::vector<double>& materials) = 0;
protected:
	std::vector<double> sigma;

	size_t row, col;
};

