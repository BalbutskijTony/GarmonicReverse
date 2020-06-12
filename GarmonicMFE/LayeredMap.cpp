#include "LayeredMap.h"


void LayeredMap::loadFromFile(const std::string& filename) {
	std::ifstream in(filename);
	double tmp;
	while (!in.eof()) {
		in >> tmp;
		sigma.push_back(tmp);
	}
	in.close();
}

// Row dont use
void LayeredMap::setRowCount(size_t count) {
	row = count;
}

void LayeredMap::setColCount(size_t count) {
	col = count;
}

double LayeredMap::getSigma(size_t index) const {
	return sigma[index / col];
}

double LayeredMap::getSigma(size_t row, size_t col) const {
	return sigma[row];
}

void LayeredMap::setSigma(size_t index, double sigma) {
	this->sigma[index / col] = sigma;
}

void LayeredMap::setSigma(size_t row, size_t col, double sigma) {
	this->sigma[row] = sigma;
}

void LayeredMap::restructMaterial(const std::vector<double>& materials) {
	for (size_t i = 0; i < sigma.size(); i++) {
		sigma[i] = materials[i];
	}
}