use strict;

my %UniqStat;
my %CigarStat;
my @InsertData;
my $Total_paired_reads = 0;
my $Total_reads_num = 0;
my $Total_bases_num = 0;

my $All_mapped_bases = 0;
my $Uniq_mapped_bases = 0;


while (<>) {
	next if(/^\@/);
	$Total_reads_num ++;
	##print $_;
	my @t = split /\t/;
	my $tag = $t[11];
	my $cigar = $t[5];
	my $equal = $t[6];
	my $insert = $t[8];
	$Total_bases_num += length($t[9]);
	
	##use the XT:A:X tag to determine unique or repetitive mapped reads
	my $uniq;
	if (/XT:A:(\w)/) {
		$uniq = $1;
	}else{
		$uniq = '*';
	}
	$UniqStat{$uniq} ++;
	
	if ($cigar =~ /\w/) {
		while ($cigar =~ /(\d+)([^\d]+)/g) {		
			$CigarStat{$2} += $1;
			if ($2 eq "M" && $uniq eq "U") {
				$Uniq_mapped_bases += $1;
			}
			
		}
	}else{
		$CigarStat{$cigar} ++;
	}

	if ($equal eq "=" && $insert =~ /\d/ && $insert != 0) {
		$Total_paired_reads ++;
		push @InsertData, abs($insert) if($Total_paired_reads < 20000 && abs($insert) < 1000);
	}

}

my $paired_reads_ratio = $Total_paired_reads / $Total_reads_num;
print "Total reads:   $Total_reads_num\n";
print "Total bases:   $Total_bases_num\n\n";

print "Statistis on reads level:\n";
foreach my $str (sort keys %UniqStat) {
	my $num = $UniqStat{$str};
	my $ratio = $num / $Total_reads_num;
	print "Unique status $str: $num  $ratio\n";
}
print "Paired reads:  $Total_paired_reads  $paired_reads_ratio\n\n";

print "Statistis on bases level:\n";
foreach my $str (sort keys %CigarStat) {
	my $num = $CigarStat{$str};
	my $ratio = $num / $Total_bases_num;
	print "Cigar status $str: $num  $ratio\n";
}

my $ratio = $Uniq_mapped_bases / $Total_bases_num;
print "\nUniq Mapped bases: $Uniq_mapped_bases  $ratio\n";

print "\nInsert size statistics:\n";
my ($total_insert, $min_insert, $max_insert, $mean_insert, $median_insert, $sd_insert);
number_stat(\@InsertData,\$total_insert, \$min_insert, \$max_insert, \$mean_insert, \$median_insert, \$sd_insert);
print "Mean:   $mean_insert\n";
print "Median: $median_insert\n";
print "Min:    $min_insert\n";
print "Max:    $max_insert\n";
print "Sd:     $sd_insert\n";


sub number_stat{
        my ($nums_p, $total_p, $min_p, $max_p, $mean_p, $median_p, $SD_p) = @_;
        
        my $add_square = 0;
        @$nums_p=sort {$a<=>$b} @$nums_p;
        my $loop=0;
        foreach  (@$nums_p) {
                $$total_p+=$_;
                $add_square+=$_*$_;
                $loop++;
        }
        my $number = @$nums_p;
        $$min_p = $nums_p->[0];
        $$max_p = $nums_p->[$number-1];
        $$mean_p = $$total_p/$number;
        $$median_p = $nums_p->[int $number/2];
        $$SD_p=sqrt( ($add_square-$$total_p*$$total_p/$number)/ ($number-1) ) if($number>1);
}

