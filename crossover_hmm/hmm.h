#ifndef __SEQ_KMER_H_
#define __SEQ_KMER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;

//load the HMM parameters from file into memory
void Load_hmm_para(string &hmm_file, int &k, int &m, vector <string> &state, vector <string> &alpha, vector <double> &initial, vector <double> &trans, vector <double> &emit);

//infer the hidden sequence by viterbi dynamic programming
void Viterbi_DP(string &observe_seq, string &hidden_seq, vector <string> &state, vector <string> &alpha, vector <double> &initial, vector <double> &trans, vector <double> &emit);


#endif
