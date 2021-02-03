#!/usr/bin/perl

# This program takes a list of SNPs and calculates the heterozygosity in windows of chosen size

# open input and output files
open (INPUTSNPFILE, "snp_locations_AD_new.txt");
open (OUTPUTFILE, "> heterozygosity_10kb_windows_new.txt");

# This defines the size of windows we will look at hets for
my $windowsize=10000;
my $halfwindowsize=$windowsize/2;

my $numchroms=14;

# ignore first line
my $header = <INPUTSNPFILE>;

$counter=0;


# Hash of chromosome lengths
$chrlength{1} = 2291499;
$chrlength{2} = 1621675;
$chrlength{3} = 1575141;
$chrlength{4} = 1084805;
$chrlength{5} = 1814975;
$chrlength{6} = 1422463;
$chrlength{7} = 1399503;
$chrlength{8} = 1398693;
$chrlength{9} = 1186808;
$chrlength{10} = 1059964;
$chrlength{11} = 1561994;
$chrlength{12} = 774062;
$chrlength{13} = 756744;
$chrlength{14} = 926563;

# read in list of snps in columns, number and cumulative location
while (<INPUTSNPFILE>)
{
        chomp($_);
		@line = split(/\t/, $_);
# We have 4 columns: Chromosome, SNP number, Position, Cum_Position in that order
		for ($a=0;$a<4;$a++)
		{
			$snps[$counter][$a] = $line[$a];
		}
		$counter++;
#		print "\n$counter";
}

$maxsnps=$counter;

print OUTPUTFILE "Chrom\tMidpoint\tHet\n";

for ($a=1;$a<=$numchroms;$a++)
{
	$end = $windowsize;
	$mid = $halfwindowsize;
	$begin = 0;
	print "\nmid=$mid chrlength($a)=$chrlength{$a}";
	until ($mid > $chrlength{$a})
	{	
		$sum=0;
#		print "\n got in";
		for ($b=0;$b<$maxsnps;$b++)
		{
			if (($snps[$b][0] == $a) && ($snps[$b][1] <= $end) && ($snps[$b][1] > $begin))
			{
				$sum++;
#				print "inside summing ";
			}
		}
		$het = $sum/$windowsize;
		print "$a\t$mid\t$het\n";
		print OUTPUTFILE "$a\t$mid\t$het\n";		
		$end += $windowsize;
		$mid += $windowsize;
		$begin += $windowsize;
	}

}


			
		


