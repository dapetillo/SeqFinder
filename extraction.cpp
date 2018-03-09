#include <algorithm>
#include "extraction.h"


std::vector <std::vector<char>> ExtractSeqs::read_seqs (std::string path){
	/*it access to the folder containing the genetics sequences.
	WARNING: only *.fasta files allowed
	*/

	DIR *dirp;
	struct dirent *drnt;

	std::vector <std::vector<char>> seqs;			// If first string is "ATCCT", final result is: 
	std::vector <std::string> dirlist;				// vector = {{'A','T','C','C','T'}, {}, {}, ...} 
	std::string temp_seq;
	std::string line;


	if ((dirp = opendir(path.c_str())) != NULL){

		while ((drnt = readdir(dirp)) != NULL){

			if (drnt->d_name[0] != '.'){

				dirlist.push_back(drnt->d_name);
			}

			else{

				continue;
			}
		}

		closedir(dirp);
	}

	else {

		perror ("Wrong path given!");
	}
	
	std::sort(dirlist.begin(), dirlist.end());
	seqs.resize(dirlist.size());					// one must always allocate memory to vectors unless using push_back

	for (unsigned i=0; i<dirlist.size(); i++){

		std::ifstream seq_file (path+"/"+dirlist[i]);
		std::cout << dirlist[i] << std::endl;
		if (seq_file.is_open()){

			std::getline(seq_file, line);			// discard the ">" line

			while (std::getline(seq_file, line)){

				if (line.back() == '\r'){			// WARNING: depending on where FASTA files come from, 
													// they may contain carriage return causing wrong reading
													// of the sequence
					line.pop_back();
				}

				temp_seq += line;
			}

			if (islower(temp_seq[0])){				// some fasta files have lowercase letters: uppercase is the format
													// to follow
				transform(temp_seq.begin(), temp_seq.end(), temp_seq.begin(), ::toupper);
			}

			
			std::copy(temp_seq.begin(), temp_seq.end(), std::back_inserter(seqs[i]));
			temp_seq.clear();
			seq_file.close();
		}
	}

	return seqs;	
}	