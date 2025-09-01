#pragma once
#include "Header.h"


class Int_point
{
public:
	std::unordered_map<string, double> parameters; // все параметры
	std::array<double, 3> center;

	Int_point(const double& a, const double& b, const double& c);
};

