#include <iostream>
#include <algorithm>
#include <numeric>
#include "freqfeat.h"



std::vector <std::string>  FreqFeat::all_words (std::vector <std::string> alphabet, int k){
	/* It calculates the permutations of all the possible
	words given k.*/

	std::map <int, std::string> alphanum;									
	std::vector <std::vector <int> > permk(1);			// 1 is necessary to initialize the vector.	
	std::vector <std::string> all_w;										
	for (unsigned i=0; i<alphabet.size(); i++){

		alphanum[i] = alphabet[i];						// fill map (key: number; value: genetic letter)			
	}

	for (int i=0; i<k; i++){												

		permk[0].push_back(0);
	}	
	
	for (unsigned i=1; i<std::pow(4.0, k); i++){		//# of combinations of words

		permk.insert(permk.end(), permk[i-1]);
		permk[i].at(k-1) += 1;
		
		for (int carry=k-1; carry>0; carry--){

			if (permk[i].at(carry)>alphabet.size()-1){

				permk[i].at(carry) = 0;										
				permk[i].at(carry-1) += 1;
			}
		}
	}
	
	std::string str;
	for (unsigned i=0; i<permk.size(); i++){			//map the numbers into std::string

		for (unsigned j=0; j<permk[i].size(); j++){
			
			for (auto it=alphanum.begin(); it!=alphanum.end(); it++){
				if (permk[i].at(j)==it->first){
					str += it->second;
				}
			}
		}
		all_w.push_back(str);
		str.clear();
	}
	return all_w;
}



std::vector <std::vector<int>> FreqFeat::words_overlay (std::vector <std::vector<char>> seqs, std::vector <std::string> all_w, int k){
	/*It performs the words extraction from the genetic
	sequences and calculates their occurrences.
	*/

	std::vector <std::map<std::string, int>> ffp(seqs.size());
	std::vector <std::vector<int>> occurrences(seqs.size());

	/*Fills maps with keys as ordered in all_w and initialises their values to 0 */
	
	for (unsigned i=0; i<ffp.size(); i++){

		for (unsigned j=0; j<all_w.size(); j++){
			ffp[i][all_w[j]];
		}
	}

	/*The effective overlay analysis*/

	for (unsigned i=0; i<seqs.size(); i++){										//sequences' index
		for(unsigned j=0; j<seqs[i].size() - k + 1; j++){						//index over the words within a sequence
			std::string word(seqs[i].begin()+j, seqs[i].begin()+j+k);
			if (word.find('N') != std::string::npos)
				continue;
			else
				ffp[i][word] += 1;
		}
		
		for (auto it=ffp[i].begin(); it!=ffp[i].end(); it++){				//inserts the (sorted) occurrences in a vector

					occurrences[i].push_back(it->second);
		}
	}

	return occurrences;
}