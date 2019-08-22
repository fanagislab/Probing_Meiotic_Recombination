cohmm according to donor's haplotype and each sperm/egg SNP, infer crossover positions by Hidden Markov Model.


Here provides two HMM parameter files, one contains two states father and mother, the other contains three states father and mother and heterozygous. The two state HMM is suitable for haploid cells (sperms, polar body 2 and female pronucleus), while the three state HMM is suitable for diploid cells( diploid sperms, i.e. two sperms merged together, or polar body 1).

The recombination rate and SNP density is different for various species, you just need to modify the parameters in the HMM parameter file, and do not need to change cohmm.


./test_data/ 

cohmm input file format£º(S0110A.chr1.cover)
#Chr    Pos     FaHap   MoHap   SpGt,Score      CoverReads[Allele1,num;Allele2,num]
chr1    254186  T       C       C,13    T,0;C,1
chr1    589084  A       G       G,9     A,0;G,1
chr1    601141  T       C       C,100   T,0;C,4
chr1    602645  T       C       C,11    T,2;C,2
chr1    729679  C       G       C,26    C,1;G,0

#Chr  chromosome ID    
Pos   chromosome position
FaHap   SNP allele on the father's haplotype
MoHap   SNP allele on the mother's haplotype
SpGt,Score   Sperm/Egg SNP allel, phred scale score.
CoverReads[Allele1,num;Allele2,num],  Sequenced bases on this sperm position, and number, [not used]


cohmm has two output files£º

S0110A.chr1.cover.crossover   crossover predictions on each sperm/egg

S0110A.chr1.cover.hiddenseq   the HMM hidden states on each postions of each sperm/egg


In the testing data, S0110A£¬S0110B£¬S0110C£¬are polar body 1, polar body 2 and female pronucleus from the same oocyte.


Commands£º
../cohmm -s 20 -n 20 ../crossover_2states.hmm S0110B.chr1.cover
