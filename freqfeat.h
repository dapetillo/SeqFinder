#ifndef FREQFEAT_H
#define FREQFEAT_H

//included dependencies

#include <vector>
#include <string>
#include <map>

class FreqFeat {
	//It defines the possible words given from alphabet and k. By default, k = 4.
	//It also defines vectors where the words' occurrences are stored.
	public:

		std::vector <std::string> all_words (std::vector <std::string> alphabet = {"A", "T", "C", "G"}, int k = 4);
		std::vector <std::vector<int>> words_overlay(std::vector <std::vector<char>> seqs, std::vector <std::string> all_w,  int k = 4);
};



#endif