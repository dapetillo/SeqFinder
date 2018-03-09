#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include "extraction.h"
#include "freqfeat.h"
#include "corrcoeff.h"
#include "kendall.h"

/*The main function settles the "SeqFinder" experiment,
that is an analysis that is supposed to find a specific
sequence (probe) in another one (usually a genome),
making use of previous calculations of correlation 
coefficients and/or confidence intervals.
A for loop runs over the genome - nucleotide by 
nucleotide - extracting each time a subsequence
of the same length of the probe; then, probe and 
subsequence are correlated. If the correlation coefficient
is equal or higher than the experimental one, or is
within the CI, the position of the subsequence in the
genome is recorded; otherwise, the loop continues.

The script assumes that the first element in ss is the
probe sequence: therefore, the user should name the files
accordingly.


Variables
---------
ss: it contains the probe and the genome sequences.

all_w: permutations of the possible words given k
	   (by default, k = 4)

vecs: copy of the probe for further analysis
counting_probe: file to record positions of probe
				in genome

subseq: string where the subsequence of the genome
		is recorded

inf_limit: experimental correlation coefficient

*/
int main(){

	ExtractSeqs ext; 
	FreqFeat hist; 
	//CorrCoeff corrs;
	Kendall mer;

	//ExtractSeqs
	std::string path = "/path/to/sequences";
	std::vector <std::vector <char>> ss = ext.read_seqs(path);

	
	//FreqFeat
	std::vector <std::string> all_w = hist.all_words();
	std::vector <std::vector<char>> vecs;
	vecs.push_back(ss[0]);
	std::ofstream counting_probe ("Counting_seq_to_find.txt");
	if (counting_probe.is_open()){

		counting_probe << "#Initial position of probe in seq (0 indexing)\n";
	}

	for (unsigned i=0; i<ss[1].size() - ss[0].size() + 1; i++){
		std::string subseq(ss[1].begin() + i, ss[1].begin() + i + ss[0].size());
		std::vector<char> temp_v(subseq.begin(), subseq.end());
		vecs.push_back(temp_v);
		std::vector <std::vector<int>> ffp = hist.words_overlay(vecs, all_w);

		//CorrCoeff	
		//std::pair <std::vector<double>, std::vector <double>> ranks = corrs.ranking(ffp[0], ffp[1]);
		//double spearman = corrs.corrcoeff(ranks.first, ranks.second);
		//double f_c = corrs.ffifc(ffp[0], ffp[1]);
		long double tau = mer.tau_b(ffp[0], ffp[1]);
		long double inf_limit = 0.648689;
		if (counting_probe.is_open() && tau >= inf_limit){

			counting_probe << int(i) << std::endl;
		} 
		vecs.erase(vecs.begin() + 1);
		temp_v.clear();
	}
	counting_probe.close();	
	
	return 0;
}
