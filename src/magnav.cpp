/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-07-21 14:38
#
# Filename:		magnav.cpp
#
# Description:
#
************************************************/

#include <iostream>
#include <unordered_set>

#include "mag_compensation/mag_compensation.hpp"

using namespace magnav;

	enum TLterm
	{
		PERMANENT,
		INDUCED,
		EDDY,
		FDM,
		BIAS
	};

void print(const auto& set)
{
    for (const auto& elem : set)
        std::cout << elem << ' ';
    std::cout << '\n';
}

int main(int argc, char *argv[]) {

	std::unordered_set<TLterm> terms={PERMANENT,INDUCED,EDDY};
	print(terms);

	std::vector<double> B(10);
	print(B);

  return 0;
}
