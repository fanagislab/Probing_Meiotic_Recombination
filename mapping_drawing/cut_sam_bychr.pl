#!/usr/bin/perl 

use File::Basename qw(basename dirname); 

if (@ARGV != 5){
	die "Usage: perl bychr.pl <ref.fai> <in.sam>  <outdir> <prefix> <samtools>\n";
}

my $ref_fai = shift @ARGV;
my $sam_file = shift @ARGV;
my $outdir = shift @ARGV;
my $samtools = shift @ARGV;
my $prefix = shift @ARGV;



my $all_result_files;
mkdir($outdir);

my @chrs;
open IN, "$ref_fai" || die "fail $ref_fai";
while (<IN>) {
	my $chr = $1 if(/^(\S+)/); ##chr
	#print STDERR $chr."\n";
	push @chrs, $chr;
}
close IN;
push @chrs, "chrUn";

foreach my $out (@chrs) {
	open $out,"| $samtools view -bS -t $ref_fai - > $outdir/$prefix.$out.bam";
}

open IN, $sam_file or die $!;

while (<IN>){
	if (/^@/){
		if (/SN:(\S+)/){
			my $out = $1;
			print $out $_;
		}
		elsif (/^\@PG/){
			foreach my $out (@chrs) {
				print $out $_;
			}
		}	
	}
	my $out = $1 if(/^\S+\s+\S+\s+(\S+)/);
	$out = "chrUn" if ($out eq "*");
	print $out $_;
}

close IN;

foreach my $out (@chrs) {
	close $out;
}

