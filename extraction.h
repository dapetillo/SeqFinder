#ifndef EXTRACTSEQS_H
#define EXTRACTSEQS_H

//included dependencies

#include <iostream>
#include <fstream>
#include <dirent.h>	
#include <vector>
#include <string>


class ExtractSeqs {
	//It defines the extraction of sequences from *.fasta file given a path.
	public:

		std::vector <std::vector <char>> read_seqs(std::string path);

};

#endif