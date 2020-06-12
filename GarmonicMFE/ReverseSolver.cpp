#include "ReverseSolver.h"

// 3
void ReverseSolver::startState() {
	if (!jacobi.empty()) jacobi.clear();
	if (!right.empty()) right.clear();

	size_t dimParameters = curParams.size();
	size_t dimEpsilon = originSolve.size();

	//curTask->init(); // ???
	curTask->setParameters(curParams);
	curTask->solve();
	curSolve = curTask->getReceiveValues();

	for (size_t i = 0; i < dimParameters; i++) {
		shiftSolve.push_back({});
		shiftParams.push_back(0.0001);//curParams[i] / 1000.);
		jacobi.push_back(curParams);
		right.push_back(0.0);
	}

	for (size_t i = 0; i < dimEpsilon; i++) {
		if (originSolve[i] != 0.0)
			weight.push_back(1. / originSolve[i]);
		else
			weight.push_back(1.0);
	}

	curIter = 0;
}

void ReverseSolver::iterate() {

}

bool ReverseSolver::isEnd() {
	return (curIter >= maxIter) || (curError <= error);
}

// 2
void ReverseSolver::setTask(Task* newTask) {
	curTask = newTask;
}


void ReverseSolver::solve() {
	curTask->setParameters(curParams);
	curTask->solve();
	curSolve = curTask->getReceiveValues();
	curError = calcError(curSolve);

	Eigen::MatrixXd Jr(curSolve.size(), curParams.size());
	Eigen::MatrixXd Jrt(curParams.size(), curSolve.size());
	Eigen::MatrixXd J_1(curParams.size(), curParams.size());
	Eigen::VectorXd uToSig(curParams.size());
	Eigen::VectorXd Ur(curSolve.size());

	// Start iterate process
	// Если она больше критической
	std::ofstream out("outReverse.txt");
	out << curIter << '\t' << curError << '\t' << curParams[0] << std::endl;
	while (!isEnd()) {
		auto tempP = curParams;
		for (size_t curP = 0; curP < curParams.size(); curP++) {
			tempP[curP] += shiftParams[curP];
			curTask->setParameters(tempP);
			curTask->solve();
			for (size_t col = 0; col < curSolve.size(); col++) {
				Jr(col, curP) = (curTask->getReceiveValues()[col] - curSolve[col]) / shiftParams[curP];
			}
			tempP[curP] -= shiftParams[curP];
		}
		Jrt = Jr.transpose();

		for (size_t i = 0; i < curSolve.size(); i++)
			Ur[i] = originSolve[i] - curSolve[i];

		J_1 = (Jrt * Jr).inverse();
		uToSig = J_1 * Jrt * Ur;
		for (size_t i = 0; i < curParams.size(); i++)
			curParams[i] -= uToSig[i];

		curTask->setParameters(curParams);
		curTask->solve();
		curSolve = curTask->getReceiveValues();
		curError = calcError(curSolve);

		curIter++;
		out << curIter << '\t' << curError << '\t' << curParams[0] << std::endl;
	}
	out.close();
}


void ReverseSolver::oldSolve() {
	// Zero iteration
	// Âûñ÷èòûâàåì òåêóùóþ îøèáêó
	curTask->setParameters(curParams);
	curTask->solve();
	curSolve = curTask->getReceiveValues();
	curError = calcError(curSolve);

	// Start iterate process
	// Åñëè îíà áîëüøå êðèòè÷åñêîé
	std::ofstream out("outReverse.txt");
	out << curIter << '\t' << curError << '\t';
	for (const auto& a : curParams)
		out << a << '\t';
	out << std::endl;
	while (!isEnd()) {
		// Âûñ÷èòûâàåì çíà÷åíèå ôóíêöèîíàëà â îïîðíûõ òî÷êàõ ïî âñåì êîîðäèíàòàì
		auto tempP = curParams;
		for (size_t curP = 0; curP < curParams.size(); curP++) {
			tempP[curP] += shiftParams[curP];
			curTask->setParameters(tempP);
			curTask->solve();
			shiftSolve[curP] = curTask->getReceiveValues();
			tempP[curP] -= shiftParams[curP];
		}
		// Ñòðîèì ïðîèçâîäíûå â âèäå âåêòîðîâ
		for (size_t curP = 0; curP < curParams.size(); curP++) {
			for (size_t curSol = 0; curSol < curSolve.size(); curSol++)
				shiftSolve[curP][curSol] = (shiftSolve[curP][curSol] - curSolve[curSol]) / shiftParams[curP];
		}
		// Ñòðîèì ìàòðèöó ÿêîáè
		for (size_t row = 0; row < curParams.size(); row++) {
			for (size_t col = row; col < curParams.size(); col++) {

				jacobi[row][col] = 0.0;
				for (size_t curSol = 0; curSol < curSolve.size(); curSol++)
					jacobi[row][col] += shiftSolve[row][curSol] * shiftSolve[col][curSol] * weight[curSol] * weight[curSol];
				if (row != col)
					jacobi[col][row] = jacobi[row][col];
			}

			// Ñòðîèì ïðàâóþ ÷àñòü
			right[row] = 0.0;
			for (size_t curSol = 0; curSol < curSolve.size(); curSol++)
				right[row] -= weight[curSol] * weight[curSol] * (curSolve[curSol] - originSolve[curSol]) * shiftSolve[row][curSol];
		}

		// Âûñ÷èòûâàåì âåêòîð äåëüòà-ïàðàìåòðîâ
		if (1 == curParams.size())
		{
			double delta = right[0] / jacobi[0][0];
			curParams[0] += delta;
		}
		else {
			Eigen::MatrixXd mat(jacobi.size(), jacobi.size());
			Eigen::VectorXd vec(jacobi.size());
			for (size_t i = 0; i < jacobi.size(); i++) {
				for (size_t j = 0; j < jacobi.size(); j++)
					mat(i, j) = jacobi[i][j];
				vec(i) = right[i];
			}
			
			Eigen::VectorXd delta = mat.lu().solve(vec);
			for (size_t row = 0; row < curParams.size(); row++)
				curParams[row] += delta[row];
		}
		// Âûñ÷èòûâàåì â íîâîé òî÷êå ôóíêöèîíàë
		curTask->setParameters(curParams);
		curTask->solve();
		curSolve = curTask->getReceiveValues();
		curError = calcError(curSolve);

		curIter++;
		out << curIter << '\t' << curError << '\t';
		for (const auto& a : curParams)
			out << a << '\t';
		out << std::endl;
	}
	out.close();
}

// 5
std::vector<double> ReverseSolver::getActualParameters() {
	return curParams;
}

// 0
void ReverseSolver::setStartParams(const std::vector<double>& params) {
	curParams = params;
}

// 1
void ReverseSolver::setTrueEpsilon(const std::vector<double>& epsilon) {
	originSolve = epsilon;
}

void ReverseSolver::setMaxIter(size_t maxIter) {
	this->maxIter = maxIter;
}

void ReverseSolver::setError(double error) {
	this->error = error;
}

double ReverseSolver::calcError(const std::vector<double>& solution) const {
	double result = 0.0;
	for (size_t i = 0; i < originSolve.size(); i++) {
		result = std::max((solution[i] - originSolve[i]) * (solution[i] - originSolve[i]), result); /* weight[i] * weight[i]*/;
	}
	return result;
}