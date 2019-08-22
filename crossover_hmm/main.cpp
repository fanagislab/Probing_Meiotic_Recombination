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

int SNP_score_cutoff = 20;
int Haplo_SNPnum_cutoff = 20;

void usage() 
{	cout << "cohumm <hmm_para.txt> <sperms_phased_snp_file>" << endl;
	cout << "   -s <int>   the score cutoff for each oocyte/sperm SNP, default=" << SNP_score_cutoff<< endl;
	cout << "   -n <int>   minimum SNP number for haplotype around crossover, default=" << Haplo_SNPnum_cutoff << endl;
	cout << "   -h          get help information" << endl;
	exit(0);
}

void Load_snp_file(string &file, string &seq, vector<string> &pos, string &chr)
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
		
		if (lineVec[4] == "-")
		{	continue;
		}
		int score = boost::lexical_cast<int>(lineVec[5]);
		if (score < SNP_score_cutoff)
		{	continue;
		}
		if (lineVec[4] == lineVec[2])
		{	seq.push_back('f');
			pos.push_back(lineVec[1]);
		}else if (lineVec[4] == lineVec[3])
		{	seq.push_back('m');
			pos.push_back(lineVec[1]);
		}
		
		chr = lineVec[0];
		
	}
}


int main(int argc, char *argv[])
{	
	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "s:n:h")) !=-1) {
		switch(c) {
			case 's': SNP_score_cutoff = atoi(optarg); break;
			case 'n': Haplo_SNPnum_cutoff = atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 3) usage();
	string hmm_para_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string sperm_phased_snp_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string sperm_hidden_file = sperm_phased_snp_file + ".hiddenseq";
	string sperm_crossover_file = sperm_phased_snp_file + ".crossover";

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

	cerr << "Alphabet(" << K << "):";
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
	
	string ObserveSeq;
	vector<string> ChrPos;
	string HiddenSeq;
	string ChrID;
	
	Load_snp_file(sperm_phased_snp_file, ObserveSeq, ChrPos, ChrID);
	Viterbi_DP(ObserveSeq, HiddenSeq, States, Alphabet, Initialization, Transition, Emission);
	
	
	//cout << ObserveSeq << endl;
	//cout << HiddenSeq << endl;

	//find crossover sites
	vector<int> CrossSites;
	vector<char> UpHaploType;
	vector<char> DownHaploType;
	CrossSites.push_back(0);
	UpHaploType.push_back(0);
	DownHaploType.push_back(0);
	for (int i=1; i<HiddenSeq.size(); i++)
	{	if (HiddenSeq[i] != HiddenSeq[i-1])
		{	//cout << ChrID << "\t" <<  ChrPos[i-1] << "\t" <<  ChrPos[i] << "\t" << HiddenSeq[i-1] << "\t" << HiddenSeq[i] << endl;
			CrossSites.push_back(i);
			UpHaploType.push_back(HiddenSeq[i-1]);
			DownHaploType.push_back(HiddenSeq[i]);
		} 
	}
	CrossSites.push_back(ObserveSeq.size()-1);
	
	ofstream crossoverfile;
	crossoverfile.open(sperm_crossover_file.c_str());
	if ( ! crossoverfile )
	{	cerr << "fail to open input file" << sperm_crossover_file << endl;
		exit(0);
	}

	crossoverfile << "#ChrID\tcrossover_pos\tcrossover_interval\tupstream_border\tdownstream_border\tcuration_status\tupHaplo_Type\tupHaplo_SNPnum\tdownHaplo_Type\tdownHaplo_SNPnum\n";
	for (int i=1; i<CrossSites.size()-1; i++)
	{	int this_idx = CrossSites[i];
		int up_idx = CrossSites[i-1];
		int down_idx = CrossSites[i+1];
		int up_haplo_SNP_num = this_idx - up_idx;
		int down_haplo_SNP_num = down_idx - this_idx ;
		int upstream_border = boost::lexical_cast<int>(ChrPos[this_idx-1]);
		int downstream_border = boost::lexical_cast<int>(ChrPos[this_idx]);
		int crossover_pos = (int)((upstream_border + downstream_border)/2);
		int crossover_interval = downstream_border - upstream_border;
		int curation_status = -1;
		if (up_haplo_SNP_num >= Haplo_SNPnum_cutoff && down_haplo_SNP_num >= Haplo_SNPnum_cutoff) 
		{	curation_status = 1;
		}	
		crossoverfile << ChrID << "\t" <<  crossover_pos <<  "\t" << crossover_interval << "\t"<< upstream_border << "\t" <<  downstream_border << "\t" << curation_status << "\t" << UpHaploType[i] << "\t" << up_haplo_SNP_num << "\t" << DownHaploType[i] << "\t" << down_haplo_SNP_num << endl;
	}
	

	//cout << ObserveSeq << endl;
	//cout << HiddenSeq << endl;
	ofstream hiddenfile;
	hiddenfile.open(sperm_hidden_file.c_str());
	if ( ! hiddenfile )
	{	cerr << "fail to open input file" << sperm_hidden_file << endl;
		exit(0);
	}
	
	hiddenfile << "ChrID\tSNPindex\tPosition\tObserveSeq\tHiddenSeq\n";
	for (int i=0; i<ObserveSeq.size(); i++)
	{	hiddenfile << ChrID << "\t" << i+1 << "\t" << ChrPos[i] << "\t"  <<  ObserveSeq[i] << "\t" << HiddenSeq[i] << endl;
	}

}	

