#include "FEGrid2.h"


size_t FEGrid2::size() const {
	return elements.size();
}

size_t FEGrid2::addElem(const FElement2& curElem) {
	elements.push_back(curElem);
	return elements.size() - 1;
}

void FEGrid2::constructFEGrid(const Grid2& geomGrid, MaterialMap* materialGrid) {
	size_t countByX = geomGrid.sizeX() - 1;
	size_t countByY = geomGrid.sizeY() - 1;

	// Можно передавать в класс экземпляр конечного элемента и строить сетку по его конструктору
	size_t curElement = 0.0;
	for (int curY = 0; curY < countByY; curY++) {
		for (int curX = 0; curX < countByX; curX++) {
			FElement2 newElem;

			newElem.setPoint(geomGrid.getNode(curX, curY), 0);
			newElem.setPoint(geomGrid.getNode(curX + 1, curY), 1);
			newElem.setPoint(geomGrid.getNode(curX, curY + 1), 2);
			newElem.setPoint(geomGrid.getNode(curX + 1, curY + 1), 3);
			newElem.setSigma(materialGrid->getSigma(curElement));

			addElem(newElem);
			curElement++;
		}
	}

	materials = materialGrid;
}

std::ostream& operator<<(std::ostream& out, const FEGrid2& grid) {
	out << "Finite Element Grid section start" << std::endl;
	for (auto elem : grid.elements) {
		out << elem;
	}
	out << "Finite Element Grid section end" << std::endl;
	return out;
}

const std::vector<FElement2>& FEGrid2::data() const {
	return elements;
}

void FEGrid2::restructMaterials(const std::vector<double>& sigmas) {
	materials->restructMaterial(sigmas);

	for(size_t i = 0; i < elements.size(); i++)
		elements[i].setSigma(materials->getSigma(i));
}
