#!/usr/bin/perl

use List::Util qw(sum);
# use warnings;

# This program takes a list of breakpoints and a file listing SNP locations and determines heterozygosity in windows flanking break points
# It then generates random simulated data sets so that the distributions can be compared.
# This version uses random genomic locations rather than random SNPs
# First	argument is a file listing number and location snps throughout genome
# The second argument is a list	of snp locations of breakpoints

# open input and output files
open (INPUTSNPFILE, $ARGV[0]);
open (INPUTBPFILE, $ARGV[1]);
open (OUTPUTRFILE, "> random_hets_500bp_new.txt");
open (OUTPUTAFILE, "> averages_hets_500bp_new.txt");

# This defines the size of windows we will calculate heterozygosity
my $windowsize=500;
my $halfwindowsize=$windowsize/2;

# Number of random data sets
my $randomizations = 1000;

my $maxbp = 18854554;

# ignore first line
my $header = <INPUTSNPFILE>;
# print $header;

# We will sum up lines in SNP file to get number of SNPs needed for randomizations
my $numSNPs = 0;

# This is the hash that stores the location data
my %snps;

# read in list of snps in columns, number and cumulative location
while (<INPUTSNPFILE>)
{
        chomp($_);
		@line = split(/\t/, $_);
		$snps{$line[0]} = $line[1];
		$numSNPs++;
}

# read in each line of list of breakpoints and calculate het
# put that value in array @hets
my @hets = ();
my @simhets = ();

# We will sum up lines in BP file to get number of breakpoints needed for randomizations
$numBPs = 0;

while (<INPUTBPFILE>)
{
    chomp($_);
    $snp=$_;
	@fr_hets = gethets($snp);
	push (@hets, @fr_hets);
	$numBPs++;
}
$mean_real = mean(@hets);
print OUTPUTAFILE "Observed mean heterozygosity: $mean_real\n";


# Now generate the heterozygosity of randomly chosen SNPs

for ($a=0; $a<$randomizations; $a++)
{
	my @randomhets = ();
	for ($b=0; $b<$numBPs; $b++)
	{
# choose random SNP
		my $random_bp = int(rand($maxbp))+1;
		@fr_r_hets = getrandhets($random_bp);
		push (@randomhets, @fr_r_hets);
	}
	$mean_random = mean(@randomhets);
	push (@simhets, $mean_random);
	print "mean random $a = $mean_random\n";
	print OUTPUTRFILE "$mean_random\n";
}

$mean_sim = mean(@simhets);
print OUTPUTAFILE "Simulated mean heterozygosity: $mean_sim\n";

sub gethets {
	my @localhet = ();
	my $focalSNP = shift;
# Move back from the SNP until we get to a distance larger than $halfwindowsize
	$end = 'no';
	$sumhets = 0;
	$leftsnp = $focalSNP-1;
	if ($snps{$focalSNP} < $halfwindowsize) {$edge="T";$denominator=$snps{$focalSNP}} else {$edge="F";}
	until ($end eq 'yes')
	{
		if ((($snps{$focalSNP} - $snps{$leftsnp} < $halfwindowsize)) && ($snps{$leftsnp}))
		{
			$sumhets++;
			$leftsnp--;
		}
		else
		{
			$end = 'yes';
		}
	}
	if ($edge eq "T") {if ($denominator != 0) {$rate = $sumhets/$denominator;} else {$rate=0;}}
	else
	{$rate = $sumhets/$halfwindowsize;}
	push (@localhet, $rate);
# Now move forward from the SNP until we get to a distance larger than $halfwindowsize
	$end = 'no';
	$sumhets = 0;
	$rightsnp = $focalSNP+1;
	if ($snps{$numSNPs}-$snps{$focalSNP} < $halfwindowsize) {$edge="T";$denominator=$snps{$numSNPs}-$snps{$focalSNP}} else {$edge="F";}
	until ($end eq 'yes')
	{
		if ((($snps{$rightsnp} - $snps{$focalSNP} < $halfwindowsize)) && ($snps{$rightsnp}))
		{
			$sumhets++;
			$rightsnp++;
		}
		else
		{
			$end = 'yes';
		}
	}
	if ($edge eq "T") {if ($denominator != 0) {$rate = $sumhets/$denominator;} else {$rate=0;}}
	else
	{$rate = $sumhets/$halfwindowsize;}
	push (@localhet, $rate);
	return @localhet;
}

sub getrandhets {
	my @localhet = ();
	my $focalbp = shift;
# Move back from the SNP until we get to a distance larger than $halfwindowsize
	$leftend = $focalbp-$halfwindowsize;
	$sumhets = 0;
	for ($c=0;$c<$numSNPs;$c++)
	{
		if (($snps{$c} > $leftend)&&($snps{$c} <= $focalbp))
		{
			$sumhets++;
		}
	}
	$rate = $sumhets/$halfwindowsize;
	push (@localhet, $rate);
# Now move forward from the SNP until we get to a distance larger than $halfwindowsize
	$sumhets = 0;
	$rightend = $focalbp+$halfwindowsize;
	for ($c=0;$c<$numSNPs;$c++)
	{
		if (($snps{$c} >= $focalbp)&&($snps{$c} < $rightend))
		{
			$sumhets++;
		}
	}
	$rate = $sumhets/$halfwindowsize;
	push (@localhet, $rate);
	return @localhet;
}

sub mean {
    return sum(@_)/@_;
}