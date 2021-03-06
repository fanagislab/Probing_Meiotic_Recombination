1. diploidsnp  <*.mpileup>
   
   If the host of the sperms or oocytes are sequenced, then we can firstly determine the SNPs of the host, which is a diploid genome, by mapping the sequencing data to the reference genome. For a given position and bi-allele hypothesis, such as allele A and G, there will be 3 genotypes AA, GG, and AG. diploidsnp was developed to determine the genotypes in the whole genome for a dipoid individual. It takes all the mpileup files generated by samtools as input.


2. spermsnp <diploid_hetsnp_file> <sperm_pileup_file>
   
   If the SNP genotypes of the host were known, we can use them to faciliate genotype the gamete cells (sperms in male,  female pronucleus and polar body 1 and polar body 2 in female). Only the heterozygous site of the host were used to genotype the gamete cells. For example, if the host genotype is AG, then a haploid gamete genotype can only be A or G with equal probability. This is true for sperms, female pronucleus and polar body 2. However, polar body 1 is a dipoid cell, with some genomic regions being heterozygous and others homozygous. spermsnp can be applied to all these cell types and just assume them as haploid cell. Although the results for polar body 1 do not make so much senese, it can still be used to infer the crossover positions in downstream analysis. The first input file is in soapsnp format, and the second input file is in samtools pileup format.  Note that if the host heterozygous SNP are not available, population SNP from dbSNP database (rpt format file) can also be used as the first input file.
   
   
   
