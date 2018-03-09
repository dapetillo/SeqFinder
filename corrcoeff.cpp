#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include "corrcoeff.h"


/*CorrCoeff calculates the Pearson and Spearman correlations for two random variables Y=F(X),
where F() is an unknown relation. Pearson is the covariance of X and Y divided by their
standard deviations: cov(X,Y)/std(X)*std(Y).
Spearman is Pearson but applied to the ranks of X and Y.
*/
std::pair <std::vector<double>, std::vector <double>> CorrCoeff::ranking (std::vector <int> x, std::vector <int> y) {
	/*it calculates the ranks of x and y, rg(x) and rg(y) respectively. 
	rg(x) are simply the indexes of x in ascending order.
	rg(y) are calculated in two steps:
	1) Assign ranks to y as done for x;
	2) Associate each y to the corresponding x and so for rg(y) with rg(x).
	In the end the function returns the pair (rg(x), rg(y)).
	*/

	/***Ranking of x***/

	std::vector <int> x_copy = x;
	std::vector <int>::iterator itx = x_copy.begin();
	std::vector <double> rank_x;

	/* Calculation of the (averaged) ranks for x*/
	std::sort(x_copy.begin(), x_copy.end());

	while (itx!=x_copy.end()){

		int conteggio = std::count(itx, x_copy.end(), *itx);
		double average = ((itx-x_copy.begin())*2.+conteggio+1.)/2.; // if x contains N same values, an average 
		for (int i=0; i<conteggio; i++){							// of the ranks is made
			rank_x.push_back(average);
		}

		if (conteggio!=1){

			itx += conteggio;
		}
		else{

			itx++;
		}
	}

	/***Ranking of y***/

	std::vector <int> y_copy = y;
	std::vector <int>::iterator ity =y_copy.begin();
	std::vector <double> temp_rank_y;

	/*Calculation of the (averaged) ranks for y ordered in ascending order*/
	std::sort (y_copy.begin(), y_copy.end());

	while (ity!=y_copy.end()){

		int conteggio = std::count(ity, y_copy.end(), *ity);
		double average = ((ity-y_copy.begin())*2.+conteggio+1.)/2.; 	// if y contains N same values, an average 
		for (int i=0; i<conteggio; i++){								// of the ranks is made
			temp_rank_y.push_back(average);
		}

		if (conteggio!=1){

			ity += conteggio;
		}
		else{

			ity++;
		}
	}

	//Sorting of the y vector by the ascending order of x. The resulting vector is stored in occs.
	std::vector <std::pair<int,int>> occs(y.size());
	std::vector <double> rank_y;
	std::multimap <double,double> y_map; // first (key): occurrences. second (value): ranking. It allows
										 // non-unique keys!


	for (unsigned i=0; i<occs.size(); i++){
		occs[i].first = x[i];
		occs[i].second = y[i];

	}	

	std::sort(occs.begin(), occs.end(), [](const std::pair<int, int> &left, const std::pair<int, int> &right){
		return left.first < right.first;
	});

	
	for (unsigned i=0; i<y_copy.size(); i++){
		y_map.insert(std::pair <int,double> (y_copy[i], temp_rank_y[i]));
	}


	
	// The loop assigns to each element in occs the right rank and stores it in rank_y.
	for (unsigned i=0; i<occs.size(); i++){

		auto it = y_map.find(occs[i].second);
		rank_y.push_back(it->second);
	}

	return std::pair <std::vector<double>, std::vector <double>> (rank_x, rank_y);
}


/*if xvalues, yvalues are occurrences, it returns Pearson; if ranks, Spearman*/
double CorrCoeff::corrcoeff (std::vector <double> xvalues, std::vector <double> yvalues){

	double mean_x;
	double mean_y;
	double cov_xy = 0;
	double var_x = 0;
	double var_y = 0;


	mean_x = std::accumulate(xvalues.begin(), xvalues.end(), 0.0) / xvalues.size();
	mean_y = std::accumulate(yvalues.begin(), yvalues.end(), 0.0) / yvalues.size();


	for (unsigned i=0; i<xvalues.size(); i++){			// x and y values must have the same size.

		cov_xy += double((xvalues[i] - mean_x)*(yvalues[i] - mean_y));
		var_x += std::pow(double(xvalues[i] - mean_x), 2);
		var_y += std::pow(double(yvalues[i] - mean_y), 2);  
	}

	double sigma_x = std::sqrt(var_x);
	double sigma_y = std::sqrt(var_y);
	double cc = cov_xy/(sigma_x*sigma_y);

	return cc;
}


double CorrCoeff::ffifc (std::vector <int> xvalues, std::vector <int> yvalues){

	double ff;
	double div_sum = 0;

	for (unsigned i=0; i<xvalues.size(); i++){

		div_sum += double(((xvalues[i]+1.)/(yvalues[i]+1.)) - 1.);
	}

	ff = - double((1./int(xvalues.size()))*std::fabs(div_sum));

	return ff;
}