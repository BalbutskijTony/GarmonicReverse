#pragma once

#include <vector>
#include <fstream>
#include <ostream>

class Grid2 {
public:
	Grid2() = default;
	~Grid2() = default;

	void fromFile(const std::string& fileName);
	Grid2 procreate(const std::string& xKoeff, const std::string& yKoeff) const;

	size_t sizeX() const;
	size_t sizeY() const;

	std::pair<double, double> getNode(size_t xInd, size_t yInd) const;

	friend std::ostream& operator<<(std::ostream& out, const Grid2& grid);
private:
	std::vector<double> xGrid;
	std::vector<double> yGrid;

	std::vector<double> uploadGrid(const std::string& fileName, const std::vector<double>& mainGrid) const;
};
