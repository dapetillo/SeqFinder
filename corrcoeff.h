#ifndef CORRCOEFF_H
#define CORRCOEFF_H

//included dependencies

#include <utility>
#include <vector>


class CorrCoeff {
	//It defines the Spearman and Pearson correlation functions.
	public:

		std::pair <std::vector<double>, std::vector <double>> ranking (std::vector <int> x, std::vector <int> y);
		double corrcoeff (std::vector <double> xvalues, std::vector <double> yvalues);
		double ffifc (std::vector <int> xvalues, std::vector <int> yvalues);
};

#endif