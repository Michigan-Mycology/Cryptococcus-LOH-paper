#!/usr/bin/perl

# This program takes a simple snp file and finds the mode SNP type in windows to reduce the complexity for plotting

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

$numstrains=39;

$maxbase=18860000;

$windowsize=10000;

# Cols is the number of columns before the real data that are variable between simple files
$cols = 3;

$firstline = <INPUTFILE>;

$snp = 1;

print OUTPUTFILE $firstline;

while (<INPUTFILE>)
{
	$nextline = $_;
	chomp ($nextline);
	@genos = split(/\t/, $nextline);
	$genos[$numstrains+$cols] =~ s/\s+$//;
	$position[$snp] = $genos[2];
	for ($a=$cols;$a<$numstrains+$cols;$a++)
	{
		$dataset[$snp][$a] = $genos[$a];
	}
	$snp++;
}

$numsnps=$snp;

$leftbase=0;
$rightbase=$leftbase+$windowsize;
$printbase=$rightbase/2;
$bases =  0;

# Initialize values at zero
for ($b=$cols;$b<$numstrains+$cols;$b++)
{	
	$numBs{$b}= 0;$numAs{$b}= 0;$numDs{$b}= 0;$numNs{$b}= 0;
}

until ($rightbase == $maxbase)
{
	for ($a=1;$a<$numsnps;$a++)
	{
		if (($position[$a]>$leftbase) && ($position[$a]<=$rightbase))
		{
			for ($b=$cols;$b<$numstrains+$cols;$b++)
			{
				if ($dataset[$a][$b] eq 'B') {$numBs{$b}++;}
				elsif ($dataset[$a][$b] eq 'A') {$numAs{$b}++;}
				elsif ($dataset[$a][$b] eq 'D') {$numDs{$b}++;}
				else {$numNs{$b}++;}
			}
			$bases++;
		}
	}
	print OUTPUTFILE "$printbase";
	for ($a=$cols; $a<$numstrains+$cols; $a++)
	{
		if (($bases == 0) | ($bases-$numNs{$a} == 0)) {$basetoprint='N';}
		elsif ($numBs{$a}/($bases-$numNs{$a}) > 0.5)
		{
			$basetoprint = 'B';
		}
		elsif ($numAs{$a}/($bases-$numNs{$a}) > 0.5)
		{
			$basetoprint = 'A';
		}
		elsif ($numDs{$a}/($bases-$numNs{$a}) > 0.5)
		{
			$basetoprint = 'D';
		}
		else
		{
			$basetoprint = 'N';
		}
		print OUTPUTFILE "\t$basetoprint";
	}
	$leftbase += $windowsize;
	$rightbase += $windowsize;
	$printbase += $windowsize;
	$bases =  0;
	print OUTPUTFILE "\n";
	for ($b=$cols;$b<$numstrains+$cols;$b++)
	{	
		$numBs{$b}= 0;$numAs{$b}= 0;$numDs{$b}= 0;$numNs{$b}= 0;
	}
}

