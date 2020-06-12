#pragma once
#include <vector>

class Task abstract
{
public:
	virtual ~Task() {};
	
	virtual void init() = 0;
	virtual size_t paramCount() const = 0;
	virtual void setParameters(const std::vector<double>& params) = 0;
	virtual void solve() = 0;
	virtual const std::vector<double>& getReceiveValues() = 0;
};

