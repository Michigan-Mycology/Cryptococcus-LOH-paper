#!/usr/bin/perl

# This program takes a text SNP file from simple conversion and extracts out those SNPs
# of high quality that are heterozygous in the ancestors and 

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl extract_hq_snps_crypto.pl inputtxtfile outfile\n";
}

$numstrains=39;

# Cols is the number of columns before the real data that are variable between simple files
$cols = 3;

print OUTPUTALLELES "Chromosome\tPosition\tCum_Position\tCB-100-1\tCB-100-2\tCB-100-3\tCB-100-5\tCB-100-6\tCC-100-1\tCC-100-2\tCC-100-3\tCC-100-4\tCC-100-5\tCC-100-6\tCN-100-1\tCN-100-2\tCN-100-3\tCN-100-4\tCN-100-5\tCN-100-6\t",
"CW-100-1\tCW-100-2\tCW-100-3\tCW-100-4\tCW-100-5\tCW-100-6\tCY-100-1\tCY-100-2\tCY-100-3\tCY-100-4\tCY-100-5\tCY-100-6\tEM3",
"\tP1W10\tP2W10\tP3W10\tP4W10\tP5W10\tP6W10\tSSD719_vitro\tSSD719_vivo\tYSB121\n";

# In the initial analysis, we looked at both the in vitro and in vivo ancestor. The in vivo ancestor seemed to have lost some heterozygosity
# So in the final script we only filter on in vitro

while (<INPUTFILE>)
{
	@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
	$genos[$numstrains+$cols-1] =~ s/\s+$//;
# Now  need to determine if SSD719_vitro (39) is a het
	unless (($genos[39] =~ /AA/) | ($genos[39] =~ /CC/) | ($genos[39] =~ /GG/) | ($genos[39] =~ /TT/) | ($genos[39] =~ /NA/))
	{
# And haploid parents EM3 (32) and YSB121 (41) are not missing or heterozygous
			if ((($genos[32] =~ /AA/) | ($genos[32] =~ /CC/) | ($genos[32] =~ /GG/) | ($genos[32] =~ /TT/)) && 
			(($genos[41] =~ /AA/) | ($genos[41] =~ /CC/) | ($genos[41] =~ /GG/) | ($genos[41] =~ /TT/)))
			{
				print OUTPUTALLELES "$genos[0]\t$genos[1]\t$genos[2]\t";
				for ($a=3;$a<$cols+$numstrains-1;$a++)
				{
					print OUTPUTALLELES "$genos[$a]\t";
				}
				print OUTPUTALLELES "$genos[$cols+$numstrains-1]\n";
			}
	}
}
		  
