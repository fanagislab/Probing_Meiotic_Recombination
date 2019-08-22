//��������ݲο�������mapped��reads���������������genotype, �˴��������ھ�������(�������Ӻ�SNP��֪)��
//����bayes theory���Զ������Ӻ�SNP�Ͳ�������ֵģ��Ϊ������ʣ�����۲�������ÿһ��genotype�ĺ�����ʡ�
//�������Ӻ�SNP��ֻ������genotype���ҳ��ֵ�������ʸ�Ϊ0.5��
//��������ֵģ�ͣ���ģ��reads���õ�ģ����ͬ���ɷ����ز������ݶ��á�����һ����ά����ref_base x cycle x seq_base x quality��
//���ͣ�����ʵ����õ�����P(D|genotype), �˴�ref_base��Ϊgeotype, ��D������������(cycle,seq_base,quality)��
//������õ���P(cycle,seq_base,quality|ref_base), �����ֵȼ���P(cycle|ref_base) * P(seq_base,quality|cycle,ref_base)��
//��ΪP(cycle|ref_base)������Ϊ��һ��������������bayes��ʽ�ļ����в������ã�������ǰ������ԡ�
//�������Խ�P(seq_base,quality|cycle,ref_base)ֱ������bayes��ʽ��Ҳ���γ���ref_base*cycle��Ϊ�У�seq_base*quality��Ϊ�еĶ�άmatrix�ļ���
//���������long double (128bit������)�����߾��ȵı�ʾ���ʡ�
//coverage depth�������ֵΪ400,�Է������������ֵָ�����ֳ���long double���ޡ�
//ѡȡ��������genotype��Ϊ���ս���������������ת��Ϊphred quality, ������Ϊ100��

//�˰��ȡsoapsnp��ʽ��error matrix��

//Author: Fan Wei, fanweisz09@gmail.com
//Date: 2012/5/16


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"

using namespace std;

typedef long double rate_t;
map<long int,vector<rate_t> > PTi; //store the prior probability for each substitutions
map<long int, string > DiGeno;
rate_t *PdTi;   //store the distribution of error and quality
int Ref_Base_num; 
int Cycle_num;
int Seq_Base_num;
int Quality_num;
int Qual_add_val;
int Cycle_add_val;
int Ref_add_val;
int Ref_Cyc_base_qual;
string Error_profile_file;
string Out_prefix = "output";
int Phred_score_cutoff = 5;
int Ascii_Qual_Start = 33;

//�ɣ������� 4�����ãǣ� N
char Bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

void usage() 
{	cout << "spermsnp <diploid_hetsnp_file> <sperm_pileup_file>" << endl;
	cout << "   version 1.0" << endl;
	cout << "   -m <file>   the matrix file of error and quality distributions, default=exe_path/illumina.bloodDNA.soapsnp.matrix" << endl;
	cout << "   -q <int>    ASCII chracter standing for quality==0, default=" << Ascii_Qual_Start << endl;
	cout << "   -s <file>   the cutoff of score in phred scale, default=" << Phred_score_cutoff << endl;
	cout << "   -h          get help information" << endl;
	exit(0);
}


//�ɣ��ãǣ� N�õ��������� 4
int get_base (char refBase)
{	
	int refBaseId;
	switch (refBase)
	{
		case 'A': refBaseId = 0; break;
		case 'C': refBaseId = 1; break;
		case 'G': refBaseId = 2; break;
		case 'T': refBaseId = 3; break;
		case 'a': refBaseId = 0; break;
		case 'c': refBaseId = 1; break;
		case 'g': refBaseId = 2; break;
		case 't': refBaseId = 3; break;
		default:  refBaseId = 4;
	}
	return refBaseId;
}


//�޶�freq_rate��1��1e-9֮�䣬ͳ����������1G
void load_soapsnp_error_profile (string &matrix_file)
{
	igzstream infile;
	infile.open(matrix_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << matrix_file << endl;
		exit(0);
	}
	
	int pi = 0; //PdTi��ѭ������
	int line_field_number = 18;
	vector<string> textlines;
	
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	if (lineStr[0] == '#')
		{	continue;
		}

		vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
		if (lineVec.size() != line_field_number)
		{	continue;
		}

		Quality_num = boost::lexical_cast<int>(lineVec[0]) + 1;
		Cycle_num = boost::lexical_cast<int>(lineVec[1]) + 1;
		
		textlines.push_back(lineStr);
	}

	Ref_Base_num = 4;
	Seq_Base_num = 4;
	Ref_Cyc_base_qual = Ref_Base_num * Seq_Base_num * Quality_num * Cycle_num;
	Qual_add_val = Ref_Base_num * Seq_Base_num * Cycle_num;
	Cycle_add_val = Ref_Base_num * Seq_Base_num;
	Ref_add_val = Seq_Base_num;

	PdTi = new rate_t[Ref_Cyc_base_qual];
	for (int j=0; j<textlines.size(); j++)
	{	
		vector<string> lineVec;
		boost::split(lineVec,textlines[j], boost::is_any_of(" \t\n"), boost::token_compress_on);
		int this_qual = boost::lexical_cast<int>(lineVec[0]); 
		int this_cycl = boost::lexical_cast<int>(lineVec[1]); 

		for (int i=2; i<line_field_number; i++)
		{	int k = i - 2;
			int this_ref = k / Seq_Base_num;
			int this_base = k % Seq_Base_num;
			pi = this_qual*Qual_add_val + this_cycl*Cycle_add_val + this_ref*Ref_add_val  + this_base;
			rate_t freq_rate = boost::lexical_cast<rate_t>(lineVec[i]); 
			
			if (pi < Ref_Cyc_base_qual)
			{	PdTi[pi] = (freq_rate > 1e-9L) ? freq_rate : 1e-9L;  
			}else{
				cerr << "pi exceeds its range " << Ref_Cyc_base_qual <<  endl;
				exit(0);
			}
		}
	}
	
}

//��snp�ļ��ж�ȡ�����ض�λ��ĳһ��genotype���ֵ��������
void assign_PTi_rates (string &diploid_snp_file)
{	
	igzstream infile;
	infile.open(diploid_snp_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << diploid_snp_file << endl;
	}
	
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of("\t"), boost::token_compress_on);
		
		long int refpos = boost::lexical_cast<long int>(lineVec[1]); 
		vector<rate_t> snp_prob(4, 0.0);  //ÿ��genotype�ĸ��ʳ�ʼ��Ϊ0
		string hetero;
		hetero.push_back(lineVec[5][0]);
		hetero.push_back(lineVec[9][0]);
		int snp1_id = get_base(lineVec[5][0]);
		int snp2_id = get_base(lineVec[9][0]);
		snp_prob[snp1_id] = 0.5;
		snp_prob[snp2_id] = 0.5;

		PTi[refpos] = snp_prob;
		DiGeno[refpos] = hetero;
	}

}

//��pileup��ʽ��base_str�еõ��ɾ���bases��Ϣ
void get_clean_bases(int refBase, string &raw_base_str,vector<int> &clean_base_vec)
{	
	int pos = 0;
	int len = raw_base_str.size();
	while (pos < len)
	{	
		if (raw_base_str[pos] == '^')
		{	pos += 2;
		}
		else if (raw_base_str[pos] == '+' || raw_base_str[pos] == '-')
		{	pos += 1;
			string indel_num;
			while (pos < len)
			{	if (raw_base_str[pos] >= 48 && raw_base_str[pos] <= 57 )  //ASCII: 0-9
				{	indel_num.push_back(raw_base_str[pos]);
					pos ++;
				}
				else
				{	break;	
				}
			}
			pos += atoi(indel_num.c_str());
		}
		else 
		{	
			switch (raw_base_str[pos])
			{
				case '.': clean_base_vec.push_back(refBase); break;
				case ',': clean_base_vec.push_back(refBase); break;
				case 'A': clean_base_vec.push_back(0); break;
				case 'a': clean_base_vec.push_back(0); break;
				case 'C': clean_base_vec.push_back(1); break;
				case 'c': clean_base_vec.push_back(1); break;
				case 'G': clean_base_vec.push_back(2); break;
				case 'g': clean_base_vec.push_back(2); break;
				case 'T': clean_base_vec.push_back(3); break;
				case 't': clean_base_vec.push_back(3); break;
				case 'N': clean_base_vec.push_back(4); break;
				case 'n': clean_base_vec.push_back(4); break;
				case '*': clean_base_vec.push_back(4); break;
				case '>': clean_base_vec.push_back(4); break;
				case '<': clean_base_vec.push_back(4); break;
				//�Զ�ʡ��$���ŵ��ж�
			}
			pos += 1;
		}
		
	}

}


//��pileup��ʽ��quality_str�еõ�quality��Ϣ����ת����phred-scale, 33��ʼ
void get_quality(vector<int> &qual_vec, string &qual_str)
{	
	
	for (int i=0; i<qual_str.size(); i++)
	{	int qual_val = qual_str[i] - Ascii_Qual_Start;
		qual_val = (qual_val < Quality_num) ? qual_val : Quality_num - 1;  
		qual_vec.push_back(qual_val);
	}
}

//��pileup��ʽ��cycle_str�еõ�cycle��Ϣ
void get_cycle (vector<int> &cycl_vec, string &cycl_str)
{	
	vector<string> temp_vec;
	boost::split(temp_vec,cycl_str, boost::is_any_of(","), boost::token_compress_on);
	for (int i=0; i<temp_vec.size(); i++)
	{	//int cycle_val = boost::lexical_cast<int>(temp_vec[i]);  //runs much slow
		int cycle_val = atoi(temp_vec[i].c_str()); //runs fast
		if (cycle_val > Cycle_num)
		{	cycle_val = Cycle_num;
		}
		cycl_vec.push_back(cycle_val - 1);  //cycle 0-99����100��cycles
	}
}



//����ָ���Ĳο����м��������ÿһ��genotype�ĺ������P(Ti|D)
//coverage depth�������ֵΪ400,�Է������������ֵָ�����ֳ���long double����
void calculte_genotype_probability(int refPose, vector<rate_t> &pTiDvec, vector<int> &baseVec, vector<int> &qualVec, vector<int> &cyclVec)
{	
	rate_t pD = 0.0;
	int len = (baseVec.size() < 400) ? baseVec.size() : 400;
	for (int genotype=0; genotype<4; genotype++)
	{	rate_t pTi = PTi[refPose][genotype];  
		rate_t pDTi = 1.0;
		if (pTi != 0)
		{
			for (int j=0; j<len; j++)
			{	
				int num = qualVec[j]*Qual_add_val + cyclVec[j]*Cycle_add_val + genotype*Ref_add_val +  + baseVec[j];
				if (num >= Ref_Cyc_base_qual)
				{	cerr << "num exceeds its range " << num << endl;
					exit(0);
				}
				pDTi *= PdTi[num];
			}
		}
		pTiDvec.push_back( pTi*pDTi );
		pD += pTi*pDTi;
	}

	for (int genotype=0; genotype<4; genotype++)
	{	pTiDvec[genotype] = pTiDvec[genotype] / pD;
	}

}

//find the genotype with maximum posterior likelihood and calcualte the phred-scale score (up to 100)
void find_genotype_score(vector<rate_t> &pTiDvec, char &sample_base, int &phred_score)
{	
	int max_id = 0;
	rate_t max_likelihood = 0.0;
	for (int i=0; i<pTiDvec.size(); i++)
	{	if (pTiDvec[i] > max_likelihood)
		{	max_id = i;
			max_likelihood = pTiDvec[i];
		}
	}
	
	sample_base = Bases[max_id];
	
	rate_t error_rate = 1 - pTiDvec[max_id];
	if (error_rate < 1e-10L)
	{	error_rate = 1e-10L;
	}
	phred_score = int(-10 * log10(error_rate));
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "m:q:s:h")) !=-1) {
		switch(c) {
			case 'm': Error_profile_file=optarg; break;
			case 'q': Ascii_Qual_Start=atoi(optarg); break;
			case 's': Phred_score_cutoff=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 3) usage();
	
	if (! Error_profile_file.size())
	{	string exepath = argv[0];
		int numh = exepath.rfind("/");
		Error_profile_file = exepath.substr(0,numh+1) + "illumina.bloodDNA.soapsnp.matrix";
	}
	cerr << "\nUsed error profile matrix file: " << Error_profile_file << endl;


	string diploid_snp_file = argv[optind++]; //optind, argv[optind++]˳��ָ���option�Ĳ���
	string sperm_pileup_file = argv[optind++]; //optind, argv[optind++]˳��ָ���option�Ĳ���

	clock_t time_start, time_end;
	time_start = clock();
	
	//caluate the prior probability of each genotype
	assign_PTi_rates(diploid_snp_file);
	cerr << "Load the snp file done:\n\n";

	//load the error and quality profile matrix into memory
	load_soapsnp_error_profile(Error_profile_file);
	cerr << "Dimensions of error profile:\n";
	cerr << "\n	Ref_Base_num " << Ref_Base_num << endl;
	cerr << "	Cycle_num " << Cycle_num << endl;
	cerr << "	Seq_Base_num " << Seq_Base_num << endl;
	cerr << "	Quality_num " << Quality_num << endl;
	cerr << "	PdTi array size: " << Ref_Cyc_base_qual << endl;
	cerr << "	Qual_add_val " << Qual_add_val << endl;
	cerr << "	Cycle_add_val " << Cycle_add_val << endl;
	cerr << "	Ref_add_val " << Ref_add_val << endl;
	
	time_end = clock();
	cerr << "\nLoad the error profile done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	//parse the pileup file, and calculate the likelihood of each genotypes
	igzstream infile;
	infile.open(sperm_pileup_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << sperm_pileup_file << endl;
	}
	
	cout << "#Chr\tPos\tRef\tDiGt\tSpGt\tScore\tDepth\tBases\tQuals\tCycles\t[all mapped bases]\n"; 
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of("\t"), boost::token_compress_on);
		
		long int refPos = boost::lexical_cast<long int>(lineVec[1]); 
		
		if ( ! PTi.count( refPos ) )
		{	continue;
		}

		int refBase = get_base(lineVec[2][0]);

		vector<int> baseVec;
		get_clean_bases(refBase,lineVec[4],baseVec);
		
		vector<int> qualVec;
		get_quality(qualVec, lineVec[5]);

		vector<int> cyclVec;
		get_cycle(cyclVec, lineVec[7]);

		//check the result of parsing pileup file
		if (baseVec.size() != qualVec.size() || baseVec.size() != cyclVec.size() || baseVec.size() == 0)
		{	cerr << "Error happens at : " << lineVec[4] << "\t" << lineVec[5] << "\t" << lineVec[7] << endl << endl;
			exit(0);
		}

		vector<int> FbaseVec;
		vector<int> FqualVec;
		vector<int> FcyclVec;
		
		//filter the observed data, base not N, and qual > 5;
		for (int i=0; i<baseVec.size(); i++)
		{	if (baseVec[i] != 4 && qualVec[i] >= 5) //������ֵС��5��ʱ��P(D|G)��������ȵģ����P(G|D) = P(G), ����û�з������ã������ֻ��ӳ������ʣ���ɴ���
			{	FbaseVec.push_back(baseVec[i]);
				FqualVec.push_back(qualVec[i]);
				FcyclVec.push_back(cyclVec[i]);
			}
		}
		if (FbaseVec.size() == 0)
		{	continue;
		}
		
		vector<rate_t> pTiDvec;
		calculte_genotype_probability(refPos,pTiDvec, FbaseVec, FqualVec, FcyclVec);
		
		char sample_base;
		int phred_score;
		find_genotype_score(pTiDvec,sample_base,phred_score);
		
		if (phred_score >= Phred_score_cutoff)
		{	
			cout << lineVec[0] << "\t" << lineVec[1] << "\t" << Bases[refBase] << "\t" << DiGeno[refPos] << "\t" << sample_base << "\t" << phred_score << "\t" << FbaseVec.size() << "\t";
			for (int i=0; i<FbaseVec.size(); i++)
			{	cout << Bases[FbaseVec[i]];
			}
			cout << "\t";
			for (int i=0; i<FbaseVec.size(); i++)
			{	cout << FqualVec[i] << ",";
			}
			cout << "\t";
			for (int i=0; i<FbaseVec.size(); i++)
			{	cout << FcyclVec[i] << ",";
			}
			cout << "\t" << baseVec.size() << "\t"<< lineVec[4] << "\t" << lineVec[5] << "\t" << lineVec[7] << endl;
		}
	}

	time_end = clock();
	cerr << "\nParse the pileup file and genotype calculation done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}

