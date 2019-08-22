#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Window_size, $Verbose,$Help);
GetOptions(
	"window:i"=>\$Window_size,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Window_size ||= 1000;
die `pod2text $0` if ( $Help);

my $chr_len_file = shift;
my $chr_pileup_file = shift;

my $chr_name = $1 if($chr_pileup_file =~ /(chr\w+)\./);


my %ChrLen;
my %PureLen;

open IN, $chr_len_file || die "fail $chr_len_file";
while (<IN>) {
	chomp;
	my @t  = split /\s+/;
	$ChrLen{$t[0]} = $t[1];
	$PureLen{$t[0]} = $t[2];
}
close IN;

my $chr_len = $ChrLen{$chr_name};
my $chr_len_pure = $PureLen{$chr_name};

print STDERR "chrId: $chr_name\n";
print STDERR "chrLen: $chr_len\n";
print STDERR "chrLenPure: $chr_len_pure\n";

my @CovLength;
my @CovDepth;
for (my $i=0; $i<$chr_len; $i+=$Window_size) {
	my $j = $i/$Window_size;
	$CovLength[$j] = 0;
	$CovDepth[$j] = 0;
}

open IN, "gzip -dc $chr_pileup_file|" || die "fail $chr_pileup_file";
while (<IN>) {
	my ($pos,$depth) = ($1,$2) if(/^\S+\s+(\d+)\s+\w\s+(\d+)/);
	my $j = ($pos - $pos % $Window_size) / $Window_size;
	if ($depth >= 1) {
		$CovLength[$j] ++;
		$CovDepth[$j] += $depth;
	}
}
close IN;

##output the result
my $total_covered_len = 0;
my $total_covered_bases = 0;
for (my $i=0; $i<$chr_len; $i+=$Window_size) {
	my $start_pos = $i + 1;
	my $j = $i/$Window_size;
	$total_covered_len += $CovLength[$j];
	$total_covered_bases += $CovDepth[$j];
	print "$chr_name\t$i\t$Window_size\t$CovLength[$j]\t$CovDepth[$j]\n";
}

my $coverage_ratio = $total_covered_len / $chr_len_pure;
my $average_depth = $total_covered_bases / $chr_len_pure;


print STDERR "covered_len: $total_covered_len\n";
print STDERR "mapped_bases: $total_covered_bases\n";
print STDERR "coverage_ratio: $coverage_ratio\n";
print STDERR "average_depth:  $average_depth\n";

