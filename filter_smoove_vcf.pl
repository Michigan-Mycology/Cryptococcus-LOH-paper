#!/usr/bin/perl

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl filter_vcf.pl inputvcffile outputvcffile\n";
}

$numstrains=39;

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# Set this to have a minimum quality of a genotype to be used. The GQ score of GATK
$GQmin = 100;

while (<INPUTFILE>)
{
	if ($_ =~ /#/)
	{print OUTPUTFILE $_;}
	else
	{
		@columns = split(/\t/, $_);
		if ($columns[5] > 500)
		{	
			@hap1data = split(/\:/, $columns[38]);
			@hap2data = split(/\:/, $columns[47]);
			@ancvivodata = split(/\:/, $columns[46]);
			@ancvitrodata = split(/\:/, $columns[45]);
			$hap1geno = $hap1data[0];
			$hap2geno = $hap2data[0];
			$ancvivogeno = $ancvivodata[0];
			$ancvitrogeno = $ancvitrodata[0];
			unless (($hap1geno =~ /1/) | ($hap2geno =~ /1/) | ($ancvivogeno =~ /1/) | ($ancvitrogeno =~ /1/))
			{
				print OUTPUTFILE $_;
			}
		}
	}
}
			
