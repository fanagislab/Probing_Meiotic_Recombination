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
#include "hmm.h"

using namespace std;

//For window size 200kb, 10 windows means 2 mega-bases resoultion
int CNV_winnum_cutoff = 10;

void usage() 
{	cout << "cnvhmm <hmm_para.txt> <normalized_coverage.tab>" << endl;
	
	cout << "   -h          get help information" << endl;
	exit(0);
}

void Load_coverage_file(string &file, vector<string> &seq, vector<string> &pos, string &chr)
{
	ifstream infile;
	infile.open(file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << file << endl;
		exit(0);
	}
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	if (lineStr[0] == '#' || lineStr[0] == ' ' || lineStr[0] == '\t' || lineStr[0] == '\n')
		{	continue;
		}
		vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of(";, \t\n"), boost::token_compress_on);
		
		pos.push_back(lineVec[2]);
		seq.push_back(lineVec[3]);
		chr = lineVec[1];
		
	}
}


int main(int argc, char *argv[])
{	
	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "n:h")) !=-1) {
		switch(c) {
			case 'n': CNV_winnum_cutoff = atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 3) usage();
	string hmm_para_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string normalized_coverage_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string cnv_hidden_file = normalized_coverage_file + ".hiddenseq";
	string cnv_block_file = normalized_coverage_file + ".cnv";

	cerr << hmm_para_file << "\n";

	vector <string> States;
	vector <string> Alphabet;
	vector <double> Initialization;
	vector <double> Transition;
	vector <double> Emission;
	int K; //number of states
	int M; //number of alphabets
	Load_hmm_para(hmm_para_file, K, M, States, Alphabet, Initialization, Transition, Emission);
	
	cerr << "States(" << K << "):";
	for (int i=0; i<States.size(); i++)
	{	cerr << "  " << States[i] ;
	}
	cerr << "\n";

	cerr << "Alphabet(" << M << "):";
	for (int i=0; i<Alphabet.size(); i++)
	{	cerr << "  " << Alphabet[i] ;
	}
	cerr << "\n";

	cerr << "Initialization(" << K << "):";
	for (int i=0; i<Initialization.size(); i++)
	{	cerr << "  " << Initialization[i] ;
	}
	cerr << "\n";
	
	cerr << "Transition(" << K*K << "):";
	for (int i=0; i<Transition.size(); i++)
	{	cerr << "  " << Transition[i] ;
	}
	cerr << "\n";
	
	cerr << "Emission(" << K*M << "):";
	for (int i=0; i<Emission.size(); i++)
	{	cerr << "  " << Emission[i] ;
	}
	cerr << "\n";
	
	//偶尔作弊的casino实测试数据
	//string ObserveSeq = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";
	//string HiddenSeq;
	//string RefSeq = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLLLLLFFFFFFFFFFFFLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLLLLLLFFFFFFFFFFF";
	//Viterbi_DP(ObserveSeq, HiddenSeq, States, Alphabet, Initialization, Transition, Emission);
	
	vector<string> ObserveSeq;
	vector<string> ChrPos;
	vector<string> HiddenSeq;
	string ChrID;
	
	Load_coverage_file(normalized_coverage_file, ObserveSeq, ChrPos, ChrID);
	Viterbi_DP(ObserveSeq, HiddenSeq, States, Alphabet, Initialization, Transition, Emission);
	
	
	//cout << ObserveSeq << endl;
	//cout << HiddenSeq << endl;

	//find crossover sites
	vector<int> CrossSites;
	vector<string> UpHaploType;
	vector<string> DownHaploType;
	CrossSites.push_back(0);
	UpHaploType.push_back("");
	DownHaploType.push_back("");
	for (int i=1; i<HiddenSeq.size(); i++)
	{	if (HiddenSeq[i] != HiddenSeq[i-1])
		{	//cout << ChrID << "\t" <<  ChrPos[i-1] << "\t" <<  ChrPos[i] << "\t" << HiddenSeq[i-1] << "\t" << HiddenSeq[i] << endl;
			CrossSites.push_back(i);
			UpHaploType.push_back(HiddenSeq[i-1]);
			DownHaploType.push_back(HiddenSeq[i]);
		} 
	}
	CrossSites.push_back(ObserveSeq.size());
	UpHaploType.push_back(HiddenSeq[ObserveSeq.size()-1]);
	
	ofstream crossoverfile;
	crossoverfile.open(cnv_block_file.c_str());
	if ( ! crossoverfile )
	{	cerr << "fail to open input file" << cnv_block_file << endl;
		exit(0);
	}

	crossoverfile << "#ChrID\tState\tIdxStart\tIdxEnd\tIdxLen\tPosStart\tPosEnd\tPosLen\n";
	for (int i=1; i<CrossSites.size(); i++)
	{	int this_idx = CrossSites[i];
		int up_idx = CrossSites[i-1];
		
		int block_win_num = this_idx - up_idx;
		string block_type = UpHaploType[i];

		int up_border = boost::lexical_cast<int>(ChrPos[up_idx]);
		int down_border = boost::lexical_cast<int>(ChrPos[this_idx-1]); 
		int block_base_length = down_border - up_border;

		crossoverfile << ChrID  << "\t" << block_type << "\t" <<  up_idx <<  "\t" << this_idx-1 << "\t"<< block_win_num << "\t" << up_border << "\t" <<  down_border << "\t" << block_base_length << endl;
	}
	

	//cout << ObserveSeq << endl;
	//cout << HiddenSeq << endl;
	ofstream hiddenfile;
	hiddenfile.open(cnv_hidden_file.c_str());
	if ( ! hiddenfile )
	{	cerr << "fail to open input file" << cnv_hidden_file << endl;
		exit(0);
	}
	
	hiddenfile << "ChrID\tWinIndex\tStartPos\tObserveSeq\tHiddenSeq\n";
	for (int i=0; i<ObserveSeq.size(); i++)
	{	hiddenfile << ChrID << "\t" << i << "\t" << ChrPos[i] << "\t"  <<  ObserveSeq[i] << "\t" << HiddenSeq[i] << endl;
	}

}	

