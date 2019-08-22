1. Function overview
  This package was originally developed for probing meiotic recombination and aneuploidy of single sperm cells by whole-genome sequencing, and it can be also applied to all gametes, including polar body 1, polar body 2 and female pronucleus. Major functions include determining the heterozygous SNPs in the host genome, determining the genotypes for each gamete cell, inferring the crossover positions and aneuploidy events for each gamete cell, as well as draw the display figures. The package was designed in a modular style, with perl and C++ coding language. 


2. Installation
  Among the four modules, mapping_drawing only includes perl programs, all of snp_calling crossover_hmm and aneuploidy_hmm has a Makefile in each module subdirectory, you just need to "cd " into the subdirectory for each module and type "make", and all the excutables will be generated.


3. Module and usage

  a. In mapping_drawing, there are programs for mapping the single cell human sperm and oocyte sequencing data to the human reference genome, by filtering low quality reads, bwa/samtools, coverage depth count, as well as a set of statistic programs. 
      
      reads_mapping_human_pipeline.pl  <Sampling_1Xdata.lib>

  b. In snp_calling, diploidsnp call diploid SNPs for the host, similar to SOAPsnp; spermsnp took pre-determined hetSNPs as prior information, assigned equal prior probabilities (0.5) to the two possible alleles on each SNP site, and estimated the posterior probabilities for each SNP allele in each gamete cell by Bayesian theory.
      
      diploidsnp  <*.mpileup>
      
      spermsnp <diploid_hetsnp_file> <sperm_pileup_file>

  c. In crossover_hmm, cohmm infer crossover positions by Hidden Markov Model, according to donor's haplotype and each sperm/egg SNP. Here provides two HMM parameter files, one contains two states father and mother, the other contains three states father and mother and heterozygous. The recombination rate and SNP density is different for various species, you just need to modify the parameters in the HMM parameter file, and do not need to change cohmm.

      cohmm -s 20 -n 20 ../crossover_2states.hmm S0110B.chr1.cover
      cohmm -s 0 -n 100 ../crossover_3states.hmm S0110A.chr1.cover

  
  d. In aneuploidy_hmm, cnvhmm infers the copy number for each segments of the chromsome and identify the aneuploidy events for the gamete cells (sperms in male, polar body 1 and polar body 2 and female pronucleus in female), by Hidden Markov Model, using the normalized coverage depth  as input information.

      cnvhmm ../HMM.para  ./chr1.sort.bam.pileup-A.gz.win.depth.normalized.win200000


4. Reference 

The methods implemented by this package are well recorded in the methods and supplementary methods of two papers:

Sijia Lu*, Chenghang Zong*, Wei Fan*, Mingyu Yang*, Probing Meiotic Recombination and Aneuploidy of Single Sperm Cells by Whole Genome Sequencing. Science 338, 1627 (2012)

Yu Hou*, Wei Fan*, Liying Yan*, et al. Genome Analyses of Single Human Oocytes. Cell 155, 1492¨C1506 (2013)