#include <iostream>
#include <cmath>
#include <algorithm>
#include <utility>
#include "kendall.h"


/*The Kendall correlation function is computed using the Knight's algorithm**,
reducing its complexity time to O(n*log(n)).

**"A Computer Method for Calculating Kendall's Tau with Ungrouped Data", William
R. Knight.
*/
long int Kendall::inv = 0;	//variable that will contain the number of discordant pairs D

std::vector<int> Kendall::merge(const std::vector<int>& left, const std::vector<int>& right)  {  
	
    // Fill the resultant vector with sorted results from both vectors  
    std::vector<int> result;  
    unsigned left_it = 0, right_it = 0;

    while(left_it < left.size() && right_it < right.size())  
    {  
        // If the left value is smaller than the right it goes next  
        // into the resultant vector  
        if(left[left_it] <= right[right_it]) 
        {  
            result.push_back(left[left_it]); 
            Kendall::inv += result.size() - 1 - left_it; 			
            left_it++;  
        }  
        else          
        {  
            result.push_back(right[right_it]);  
            right_it++; 
   
        }
    }  
  
    // Pushes the remaining data from both vectors onto the resultant  
    while(left_it < left.size())  
    {  
        result.push_back(left[left_it]); 
        Kendall::inv += result.size() - 1 - left_it;
        left_it++;  
    }  
  
    while(right_it < right.size())  
    {  
        result.push_back(right[right_it]);  
        right_it++;  
    }  
  	
    return result;  
}  




std::vector<int> Kendall::mergesort(std::vector<int>& vec)  {  
    // Termination condition: List is completely sorted if it  
    // only contains a single element.  
    if(vec.size() == 1)  
    {  
        return vec;  
    }  
  
    // Determines the location of the middle element in the vector  
    std::vector<int>::iterator middle = vec.begin() + (vec.size() / 2);  
  
    std::vector<int> left(vec.begin(), middle);  
    std::vector<int> right(middle, vec.end());  
  
    // Perform a merge sort on the two smaller vectors  
    left = mergesort(left);  
    right = mergesort(right);  
  
    return merge(left, right);  
}  

//The effective computation of the Kendall coefficient.
long double Kendall::tau_b (std::vector <int> x, std::vector <int> y){

    if (x.size()!=y.size())
        return 1;

    long int pairs = x.size()*(x.size() - 1)*0.5;
    long int t_x = 0;
    long int t_y = 0;
    long int t_xy = 0;


    std::vector <std::pair <int,int>> comb(x.size());

    for (unsigned i=0; i< x.size(); i++){

        comb[i].first = x[i];
        comb[i].second = y[i];
    }

    /*it orders the y occurrences following the ascending order of x; 
    ties are broken according to y (if x_i == x_j, then y_i <= y_j). The resulting vector is stored in comb. */
    std::sort(comb.begin(), comb.end(), [](const std::pair<int, int> &left, const std::pair<int, int> &right){
        if (left.first == right.first)
            return left.second < right.second;
        else
            return left.first < right.first;
    });

    std::vector <int> x_sorted;
    std::vector <int> y_sorted_by_x;

    for (unsigned i=0; i<x.size(); i++){

        x_sorted.push_back(comb[i].first);
        y_sorted_by_x.push_back(comb[i].second);
    }
    
    std::vector <int> y_mergesort = mergesort(y_sorted_by_x);
    long int swap = Kendall::inv;

    std::vector <int>::iterator itym = y_mergesort.begin();

    /*counting of tied pairs in y*/
    while (itym != y_mergesort.end()){

        long int count = std::count(itym, y_mergesort.end(), *itym);
        t_y += (count*(count-1))/2;
        if (count!=1)
            itym += count;
        else
            itym++;
    }


    long int count_x = 1;
    long int count_xy = 1;

    /*counting of tied pairs in x and jointed pairs*/
    for (unsigned i=1; i<=x_sorted.size(); i++){
        if (x_sorted[i]==x_sorted[i-1]){
            count_x += 1;
            if (y_sorted_by_x[i]==y_sorted_by_x[i-1]){
                count_xy += 1;
            }
            else {
                t_xy += (count_xy*(count_xy-1))/2;
                count_xy = 1;
            }

        }
        else {
            t_x += (count_x*(count_x-1))/2;
            t_xy += (count_xy*(count_xy-1))/2;
            count_xy = 1;
            count_x = 1;
            
        }
    }


    /*calculation of tau_b (it reduces to tau_a if t_x = t_y = 0)*/
    long double numer = pairs - t_x - t_y + t_xy - (2*swap);
    long double denom = std::sqrt((pairs - t_x)*(pairs - t_y));
    long double tau_b = numer/denom;
    Kendall::inv = 0;

    return tau_b;
}