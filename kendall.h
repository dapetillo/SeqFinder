#ifndef KENDALL_H
#define KENDALL_H

//included dependencies

#include <vector>
#include <algorithm>
class Kendall {
	//It defines the Kendall correlation function (it takes into account pair ties).
	public:
		
		std::vector<int> mergesort(std::vector<int>& vec);
		std::vector <int> merge(const std::vector<int>& left, const std::vector<int>& right);
		long double tau_b(std::vector <int> x, std::vector <int> y);
		static long int inv;
};


#endif