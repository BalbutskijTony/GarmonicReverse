#include "Grid2.h"

std::vector<double> Grid2::uploadGrid(const std::string& fileName, const std::vector<double>& mainGrid) const {

	std::vector<double> result;

	std::ifstream in(fileName);
	for (int curElem = 0; curElem < mainGrid.size() - 1; curElem++) {
		double start = mainGrid[curElem];
		double end = mainGrid[curElem + 1];

		double n, k, h;
		in >> n >> k;
		n = n - 1;
		result.push_back(start);
		if (k != 1.) {
			h = (end - start) * (1 - k) / (1 - pow(k, n));
			for (int i = 1; i < n; i++) {
				result.push_back(result.back() + h);
				h *= k;
			}
		}
		else {
			h = (end - start) / (n);
			for (int i = 1; i < n; i++) {
				result.push_back(result.back() + h);
			}
		}
	}
	result.push_back(mainGrid.back());
	in.close();
	return result;
}

Grid2 Grid2::procreate(const std::string& xKoeff, const std::string& yKoeff) const {
	Grid2 result;
	result.xGrid = uploadGrid(xKoeff, this->xGrid);
	result.yGrid = uploadGrid(yKoeff, this->yGrid);
	return result;
}

void Grid2::fromFile(const std::string& fileName) {
	std::ifstream in(fileName);
	double tmp;
	while (!in.eof()) {
		in >> tmp;
		if(xGrid.empty() || tmp != xGrid.back())
			xGrid.push_back(tmp);
		in >> tmp;
		if (yGrid.empty() || tmp != yGrid.back())
			yGrid.push_back(tmp);
	}
	in.close();
}

size_t Grid2::sizeX() const { return xGrid.size(); }
size_t Grid2::sizeY() const { return yGrid.size(); }

std::pair<double, double> Grid2::getNode(size_t xInd, size_t yInd) const {
	return std::make_pair(xGrid[xInd], yGrid[yInd]);
}

std::ostream& operator<<(std::ostream& out, const Grid2& grid) {
	out << "Main Grid section start" << std::endl;
	out << "X:\t";
	for (auto x : grid.xGrid)
		out << x << "\t";
	out << std::endl;
	out << "Y:\t";
	for (auto y : grid.yGrid)
		out << y << "\t";
	out << std::endl;
	out << "Main Grid section end" << std::endl << std::endl;
	return out;
}