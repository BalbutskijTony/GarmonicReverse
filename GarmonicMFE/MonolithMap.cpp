#include "MonolithMap.h"

void MonolithMap::loadFromFile(const std::string& filename) {
	std::ifstream in(filename);
	double tmp;
	in >> tmp;
	sigma.push_back(tmp);
	in.close();
}

void MonolithMap::setRowCount(size_t count) {
	row = count;
}

void MonolithMap::setColCount(size_t count) {
	col = count;
}

double MonolithMap::getSigma(size_t index) const {
	return sigma[0];
}

double MonolithMap::getSigma(size_t row, size_t col) const {
	return sigma[0];
}

void MonolithMap::setSigma(size_t index, double sigma) {
	this->sigma[0] = sigma;
}

void MonolithMap::setSigma(size_t row, size_t col, double sigma) {
	this->sigma[0] = sigma;
}

void MonolithMap::restructMaterial(const std::vector<double>& materials) {
	this->sigma[0] = materials[0];
}