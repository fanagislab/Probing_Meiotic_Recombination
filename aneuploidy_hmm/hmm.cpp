#include "hmm.h"

void Load_hmm_para(string &hmm_file, int &k, int &m, vector <string> &state, vector <string> &alpha, vector <double> &initial, vector <double> &trans, vector <double> &emit)
{
	ifstream infile;
	infile.open(hmm_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << hmm_file << endl;
		exit(0);
	}
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	if (lineStr[0] == '\n' || lineStr[0] == ' ' || lineStr[0] == '\t')
		{	continue;
		}
		
		if (lineStr[0] == '#')
		{	vector<string> lineVec;
			boost::split(lineVec,lineStr, boost::is_any_of(":, \t\n"), boost::token_compress_on);
			if (lineVec[0] == "#States")
			{	for (int i =1; i<lineVec.size(); i++)
				{	state.push_back(lineVec[i]);
				}
				k = state.size();
			}

			if (lineVec[0] == "#Alphabet")
			{	for (int i =1; i<lineVec.size(); i++)
				{	alpha.push_back(lineVec[i]);
				}
				m = alpha.size();
			}
			
			if (lineVec[0] == "#Initial")
			{	for (int i =1; i<lineVec.size(); i++)
				{	
					initial.push_back( boost::lexical_cast<double>(lineVec[i]) );
				}
			}
			
			if (lineVec[0] == "#Transition")
			{	trans.resize(k*k);
				int i = 0;
				while (getline( infile, lineStr, '\n' ))
				{	
					vector<string> lineVec;
					boost::split(lineVec,lineStr, boost::is_any_of(":, \t\n"), boost::token_compress_on);
					
					if (lineStr[0] == ' ' || lineStr[0] == '\n' || lineVec[0] == "S\\S")
					{	continue;
					}
					
					for (int j = 1; j<lineVec.size(); j++)
					{	
						trans[i * k + j - 1] = boost::lexical_cast<double>(lineVec[j]);
					}
					i ++;
					if (i >= k)
					{	break;
					}
				}
			}
			
			if (lineVec[0] == "#Emission")
			{	emit.resize(k*m);
				int i = 0;
				while (getline( infile, lineStr, '\n' ))
				{
					vector<string> lineVec;
					boost::split(lineVec,lineStr, boost::is_any_of(":, \t\n"), boost::token_compress_on);
					if (lineStr[0] == ' ' || lineStr[0] == '\n' || lineVec[0] == "S\\A")
					{	continue;
					}
					
					for (int j = 1; j<lineVec.size(); j++)
					{	
						emit[i * m + j - 1] = boost::lexical_cast<double>(lineVec[j]);
					}
					i ++;
					if (i >= k)
					{	break;
					}
				}
			}
		}
		
	}
}


void Viterbi_DP(vector<string> &observe_seq, vector<string> &hidden_seq, vector <string> &state, vector <string> &alpha, vector <double> &initial, vector <double> &trans, vector <double> &emit)
{
	int K = state.size();
	int M = alpha.size();
	int J = observe_seq.size();
	
	map<string,int> AlphaOrder;
	for (int m=0; m<alpha.size(); m++)
	{	AlphaOrder[alpha[m]] = m;
		//cerr << alpha[i] << " "<< AlphaOrder[alpha[i]] << endl;
	}

	cerr << "Seq_Len: " << J << endl;
	vector<double> DPscore(K*J);
	vector<int> DPlink(K*J);
	
	//initialization
	int j = 0;
	for (int k=0; k<K; k++)
	{	DPscore[k*J+j] = log(initial[k] * emit[k*M+AlphaOrder[observe_seq[j]]]); 
		//DPscore[k*J+j] = initial[k] * emit[k*M+AlphaOrder[observe_seq[j]]]; 
		//cerr << DPscore[k*J+j] << endl;
	}

	//iterate calculating all the dynamic programming cells
	for (j=1; j<J; j++)
	{	for (int k=0; k<K; k++)
		{	double max_score;
			int max_link = -1;
			for (int p=0; p<K; p++)
			{	double score = DPscore[p*J+j-1] + log(trans[p*K+k] * emit[k*M+AlphaOrder[observe_seq[j]]]);
				//double score = DPscore[p*J+j-1] * trans[p*K+k] * emit[k*M+AlphaOrder[observe_seq[j]]];
				if (max_link == -1)
				{	max_score = score;
					max_link = p;
				}else
				{	if (score > max_score)
					{	max_score = score;
						max_link = p;
					}
				}
			}
			DPscore[k*J+j] = max_score;
			DPlink[k*J+j] = max_link;
			//cerr << j << " " << k << " " << max_score << " " << max_link << endl;
		}
	}

	//trace back and output the hidden sequence
	j = J - 1;
	int end_pos;
	double max_score;
	for (int k=0; k<K; k++)
	{	
		if (k == 0)
		{	end_pos = k;
			max_score = DPscore[k*J+j];
		}else{
			if (DPscore[k*J+j] > max_score)
			{	end_pos = k;
				max_score = DPscore[k*J+j];
			}
		}
	}
	
	while (j >= 0)
	{	hidden_seq.push_back(state[end_pos]);
		end_pos = DPlink[end_pos*J+j];
		j --;
	}
	reverse(hidden_seq.begin(), hidden_seq.end());

}
