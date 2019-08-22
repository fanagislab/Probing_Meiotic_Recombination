#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl reads_mapping_pipeline.pl <sperm_reads_files.list>
  --reference <str>   input the reference genome file, default=human,hg19.fa
  --maxInsSize <int>  max allowed insert size for PE-end , default=800
  --threadNum <int>   parallel level, number of threads or jobs, default=10
  --qsub <str>        whether use qsub mode: yes or no; default=yes
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

   nohup time perl reads_mapping_pipeline.pl testing_files.list

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help);
my ($refseqfile, $maxInsSize, $threadNum, $Qsub);
GetOptions(
	"reference:s"=>\$refseqfile,
	"maxInsSize:i"=>\$maxInsSize,
	"threadNum:i"=>\$threadNum,
	"qsub:s"=>\$Qsub,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$refseqfile ||= "/WPS/BP/fanwei/database/human_database_sperm/hg19_for_sperm.fa";
my $readsCounter = "/WPS/BP/fanwei/exe/readsCounter";
my $trim_malbac_adapter_lowqual = "/WPS/BP/fanwei/exe/trim_malbac_adapter_lowqual";
my $bwa = "/WPS/BP/fanwei/exe/bwa";
my $samtools = "/WPS/BP/fanwei/exe/samtools";
my $multi_process = "/WPS/BP/fanwei/bin/multi-process.pl ";
my $qsub_sge   = "/WPS/BP/fanwei/bin/qsub-sge.pl ";

$maxInsSize ||= 800;
$threadNum ||= 10;
$Qsub ||= "yes";
my $qsub_job_num = 100;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $reads_files_list = shift;
my $reads_files_base = basename($reads_files_list);
$reads_files_base =~ s/\.\w+$//;

my %ReadFiles;
my $sampleId;
my $SeqMode;


print STDERR "Reference file: $refseqfile\n";
print STDERR "Thread Number:  $threadNum\n";

##Load the sequencing reads files into %ReadFiles
open IN, $reads_files_list || die "fail $reads_files_list";
while (<IN>) {
	s/^\s+//;
	next if(!$_);
	if(/^SampleID:\s*(\S+)/){
		$sampleId = $1 ;
	}
	elsif(/(Pair-End|Single-End)/){
		$SeqMode = $1 ;
	}
	else{
		s/\s+$//;
		push @{$ReadFiles{$sampleId}{$SeqMode}}, $_;
	}
}
close IN;
print Dumper \%ReadFiles;

mkdir("Sequencing_Data_statistics");
my $trimming_shell = "$reads_files_base.Shell_10_All_trimming.sh";
open OUT, ">$trimming_shell" || die "fail $trimming_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	my $sample_p = $ReadFiles{$sampleId};
	foreach my $seqMode (sort keys %$sample_p) {
		my $mode_p = $sample_p->{$seqMode};
		
		for (my $i=0; $i<@$mode_p; $i++) {
			my $read_file = $mode_p->[$i];
			my $fq_file = basename($read_file);
			print OUT "$trim_malbac_adapter_lowqual  -m /WPS/BP/fanwei/exe/MALBAC_primer.fa  -i /WPS/BP/fanwei/exe/illumina_adapter.fa   -t $threadNum  $read_file  ./$fq_file.trimmed.gz  ./$fq_file.trimmed.stat; mv ./$fq_file.trimmed.stat Sequencing_Data_statistics/\n";
		}
	}
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge   --maxjob $qsub_job_num  --resource vf=1G --reqsub $trimming_shell`;
}else{
	`time sh $trimming_shell`;
}
print STDERR "\n## Finished running $trimming_shell\n\n";


##(1). bwa aligning parameters:
##-o and -e allow open one gap with max gap size 10bp
##-i do not put an indel within 15 bp towards the ends 
##-q quality threshold for read trimming down to 10bp
my $aligning_shell = "$reads_files_base.Shell_11_All_aligning.sh";
open OUT, ">$aligning_shell" || die "fail $aligning_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	my $sample_p = $ReadFiles{$sampleId};
	foreach my $seqMode (sort keys %$sample_p) {
		my $mode_p = $sample_p->{$seqMode};
		
		for (my $i=0; $i<@$mode_p; $i++) {
			my $read_file = $mode_p->[$i];
			my $fq_file = basename($read_file);
			print OUT "$bwa aln -o 1 -e 10 -i 15 -q 10  -t $threadNum $refseqfile ./$fq_file.trimmed.gz > ./$fq_file.sai;  rm ./$fq_file.trimmed.gz\n";
			print OUT "$readsCounter -o ./$fq_file.count $read_file; mv ./$fq_file.count Sequencing_Data_statistics/\n";
		}
	}
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge     --maxjob $qsub_job_num  --resource vf=10G --reqsub $aligning_shell`;
}else{
	`time sh $aligning_shell`;
}

print STDERR "\n## Finished running $aligning_shell\n\n";

my $pairing_shell = "$reads_files_base.Shell_12_All_pairing.sh";
open OUT, ">$pairing_shell" || die "fail $pairing_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	my $sample_p = $ReadFiles{$sampleId};
	foreach my $seqMode (sort keys %$sample_p) {
		my $mode_p = $sample_p->{$seqMode};

		if ($seqMode eq "Pair-End") {
			for (my $i=0; $i<@$mode_p; $i+=2) {
				my $read1file = $mode_p->[$i];
				my $fq1_file = basename($read1file);
				my $read2file = $mode_p->[$i+1];
				my $fq2_file = basename($read2file);
				print OUT "$bwa sampe -a $maxInsSize $refseqfile ./$fq1_file.sai ./$fq2_file.sai $read1file $read2file | perl $Bin/cut_sam_bychr.pl $refseqfile.fai -  ./$sampleId  $samtools  $fq1_file;  rm  ./$fq1_file.sai ./$fq2_file.sai;\n"; ##
			}
		}
		if ($seqMode eq "Single-End") {
			for (my $i=0; $i<@$mode_p; $i++) {
				my $read1file = $mode_p->[$i];
				my $fq1_file = basename($read1file);
				print OUT "$bwa samse $refseqfile ./$fq1_file.sai  $read1file | perl $Bin/cut_sam_bychr.pl $refseqfile.fai -  ./$sampleId $samtools $fq1_file; rm ./$fq1_file.sai; \n";  #
			}
		}
	}
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge     --maxjob $qsub_job_num  --resource vf=10G --reqsub $pairing_shell`;
}else{
	`time $multi_process -cpu $threadNum $pairing_shell`;
}

print STDERR "\n## Finished running $pairing_shell\n\n";



##(2) sort the individual bam files
##the unmapped reads are stored in chrUn.bam, and do not need sorting
my $sorting_shell = "$reads_files_base.Shell_13_All_sorting.sh";
my $bam_files_str;
open OUT, ">$sorting_shell" || die "fail $sorting_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	foreach my $bam_file (glob("$sampleId/*.bam")) {
		my $sort_file = $1 if($bam_file =~ /(.+)\.bam$/);
		print OUT "$samtools sort -m 200000000 $bam_file $sort_file.sort\n";
		$bam_files_str .= "  $bam_file";
	}
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge     --maxjob $qsub_job_num  --resource vf=2G --reqsub $sorting_shell`;
}else{
	`time $multi_process -cpu $threadNum $sorting_shell`;
}

`rm $bam_files_str`; ##only keep the sorted bam files
print STDERR "\n## Finished running $sorting_shell\n\n";


##(3) merge the individual bam files by chrs
my $merging_shell = "$reads_files_base.Shell_14_All_merging.sh";
open OUT, ">$merging_shell" || die "fail $merging_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	my %BamFiles;
	foreach my $bam_file (glob("$sampleId/*.sort.bam")) {
		my $chr = $1 if($bam_file =~ /([^\.]+)\.sort\.bam$/);
		push @{$BamFiles{$chr}}, $bam_file;
	}
	##当只有一个bam文件时不能用samtools merge
	foreach my $chr (sort keys %BamFiles) {
		my $chr_p = $BamFiles{$chr};
		my $input_bam_files = join(" ", @$chr_p);
		if (@$chr_p > 1) {
			print OUT "$samtools merge $sampleId/$chr.sort.bam $input_bam_files; rm $input_bam_files\n";
		}elsif(@$chr_p == 1) {
			print OUT "mv $input_bam_files $sampleId/$chr.sort.bam\n";
		}
	}
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge     --maxjob $qsub_job_num  --resource vf=2G --reqsub $merging_shell`;
}else{
	`time $multi_process -cpu $threadNum $merging_shell`;
}


##执行完此步骤，只剩下了每个sample总的按染色体分的sorted.bam文件。
print STDERR "\n## Finished running $merging_shell\n\n";


##(4) pileup and depth statistics
my @chrs;
open IN, "$refseqfile.fai" || die "fail $refseqfile.fai";
while (<IN>) {
	my $chr = $1 if(/^(\S+)/); ##chr
	push @chrs, $chr;
}
close IN;

my $pileup_shell = "$reads_files_base.Shell_15_All_pileup.sh";
open OUT, ">$pileup_shell" || die "fail $pileup_shell";
foreach my $sampleId (sort keys %ReadFiles) {
	foreach my $chrId (@chrs) {
		print OUT "$samtools view $sampleId/$chrId.sort.bam | perl $Bin/sam_mapping_stat.pl > $sampleId/$chrId.sort.bam.mapping.stat\n";
		print OUT "$samtools mpileup -d 1000000 -O -s -f $refseqfile $sampleId/$chrId.sort.bam | gzip - > $sampleId/$chrId.sort.bam.pileup.gz; perl $Bin/calcu_window_depth.pl --window 100 $refseqfile.len $sampleId/$chrId.sort.bam.pileup.gz > $sampleId/$chrId.sort.bam.pileup.gz.win.depth 2> $sampleId/$chrId.sort.bam.pileup.gz.win.depth.stat; \n";
		print OUT "$samtools mpileup -A -d 1000000 -O -s -f $refseqfile $sampleId/$chrId.sort.bam | gzip - > $sampleId/$chrId.sort.bam.pileup-A.gz; perl $Bin/calcu_window_depth.pl --window 100 $refseqfile.len $sampleId/$chrId.sort.bam.pileup-A.gz > $sampleId/$chrId.sort.bam.pileup-A.gz.win.depth 2> $sampleId/$chrId.sort.bam.pileup-A.gz.win.depth.stat; rm $sampleId/$chrId.sort.bam.pileup-A.gz\n";
	}	
}
close OUT;
if ($Qsub eq "yes") {
	`time $qsub_sge    --maxjob $qsub_job_num --resource vf=2G --reqsub $pileup_shell`;
}else{
	`time $multi_process -cpu $threadNum $pileup_shell`;
}
print STDERR "\n## Finished running $pileup_shell\n\n";


##Finally, make statistics on sequence amonunt and mapping results
my $stastistic_file = "$reads_files_base.All_samples_reads_mapping_statistics.tab";
open OUT, ">$stastistic_file" || die "fail $stastistic_file";
print OUT "sampleId\traw_reads_num\traw_bases_num\tclean_reads_num\tclean_bases_num\tclean_bases_ratio\tmapped_bases_num\tmapped_bases_ratio\tuniquely_mapped_bases_num\tuniquely_mapped_bases_ratio\tgenome_coverage_ratio\taverage_sequencing_depth\tInsert_average\tInsert_SD\tsex_type\n";
foreach my $sampleId (sort keys %ReadFiles) {
	my $raw_reads_num = 0;
	my $raw_bases_num = 0;
	my $clean_reads_num = 0;
	my $clean_bases_num = 0;
	my $mapped_bases_num = 0;
	my $mapped_bases_ratio = 0.0;
	my $coverage_ratio = 0.0;
	my $average_depth = 0.0;
	my $sample_p = $ReadFiles{$sampleId};
	foreach my $seqMode (sort keys %$sample_p) {
		my $mode_p = $sample_p->{$seqMode};
		
		for (my $i=0; $i<@$mode_p; $i++) {
			my $read_file = $mode_p->[$i];
			my $fq_file = basename($read_file);
			my $str = `cat Sequencing_Data_statistics/$fq_file.trimmed.stat`;
			$raw_reads_num += $1 if($str =~ /total_raw_reads:\s+(\d+)/);
			$raw_bases_num += $1 if($str =~ /total_raw_bases:\s+(\d+)/);
			$clean_reads_num += $1 if($str =~ /total_clean_reads:\s+(\d+)/);
			$clean_bases_num += $1 if($str =~ /total_clean_bases:\s+(\d+)/);
		}
	}

	my $str = `cat $sampleId/*pileup-A.gz.win.depth.stat`;
	my $total_pure_chrlen = 0;
	my $total_covered_len = 0;
	my $total_mapped_bases = 0;
	while ($str =~ /chrLenPure:\s+(\d+)/g) {
		$total_pure_chrlen += $1;
	}
	while ($str =~ /covered_len:\s+(\d+)/g) {
		$total_covered_len += $1;
	}
	while ($str =~ /mapped_bases:\s+(\d+)/g) {
		$total_mapped_bases += $1;
	}
	$mapped_bases_ratio = $total_mapped_bases / $clean_bases_num;
	$coverage_ratio = $total_covered_len / $total_pure_chrlen;
	$average_depth = $total_mapped_bases / $total_pure_chrlen;
	
	my $str = `cat $sampleId/chrX.sort.bam.pileup.gz.win.depth.stat`;
	my $chrX_avg_depth = $1 if($str =~ /average_depth:\s+(\S+)/);
	my $str = `cat $sampleId/chrY.sort.bam.pileup.gz.win.depth.stat`;
	my $chrY_avg_depth = $1 if($str =~ /average_depth:\s+(\S+)/);
	my $sex_type = ($chrX_avg_depth > $chrY_avg_depth) ? "X" : "Y";
	
	my $uniq_mapped_bases = 0;
	my $uniq_mapped_ratio = 0;

	my $str = `cat $sampleId/*.mapping.stat`;
	while ($str =~ /Uniq Mapped bases:\s+(\d+)/g) {
		$uniq_mapped_bases += $1;
	}
	$uniq_mapped_ratio = $uniq_mapped_bases / $clean_bases_num;
	
	my $Insert_average = 0;
	my $Insert_SD = 0;
	$Insert_average = $1 if($str =~ /Mean:\s+(\d+)/);
	$Insert_SD = $1 if($str =~ /Sd:\s+(\d+)/);
	
	my $clean_bases_ratio = $clean_bases_num / $raw_bases_num;
	print OUT "$sampleId\t$raw_reads_num\t$raw_bases_num\t$clean_reads_num\t$clean_bases_num\t$clean_bases_ratio\t$total_mapped_bases\t$mapped_bases_ratio\t$uniq_mapped_bases\t$uniq_mapped_ratio\t$coverage_ratio\t$average_depth\t$Insert_average\t$Insert_SD\t$sex_type\n";
	
	print STDERR "total_pure_chrlen: $total_pure_chrlen\n";
	print STDERR "total_covered_len: $total_covered_len\n";
	print STDERR "total_mapped_bases: $total_mapped_bases\n";
}
close OUT;

##draw coverage depth figures
my $drawing_shell = "$reads_files_base.Shell_16_All_drawing.sh";
open OUT, ">$drawing_shell" || die "fail $drawing_shell";
mkdir("Coverage_depth_figure");
foreach my $sampleId (sort keys %ReadFiles) {
	print OUT "perl $Bin/draw_sperm_coverage_haplotype_crossover.pl $sampleId ./$sampleId  ./Coverage_depth_figure  ./$reads_files_base.All_samples_reads_mapping_statistics.tab\n";
}
if ($Qsub eq "yes") {
	`time $qsub_sge     --maxjob $qsub_job_num  --resource vf=5G --reqsub $drawing_shell`;
}else{
	`time $multi_process -cpu $threadNum $drawing_shell`;
}
print STDERR "\n## Finished running $drawing_shell\n\n";


####################################################
################### Sub Routines ###################
####################################################
