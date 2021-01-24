#!/usr/bin/perl

# This program takes a vcf file and makes a file that is simplified for snp analysis

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl convert_vcf_to_simple_crypto.pl inputvcffile outfile\n";
}

$numstrains=39;

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# mindepth is the minimum number of reads to not filter out a pool
$mindepth = 18;

# maxdepthpool is the maximum number of reads or it should be filtered out
$maxdepth = 350;

# Threshold is the max number of strains for which missing data is possible for a locus to be included
$threshold = 10;

# Set this to have a minimum quality of a genotype to be used. The GQ score of GATK
$GQmin = 99;

# Hash of cumulative lengths
$chrtotal{cneoH99_Chr1} = 0;
$chrtotal{cneoH99_Chr2} = $chrtotal{cneoH99_Chr1} + 2291499;
$chrtotal{cneoH99_Chr3} = $chrtotal{cneoH99_Chr2} + 1621675;
$chrtotal{cneoH99_Chr4} = $chrtotal{cneoH99_Chr3} + 1575141;
$chrtotal{cneoH99_Chr5} = $chrtotal{cneoH99_Chr4} + 1084805;
$chrtotal{cneoH99_Chr6} = $chrtotal{cneoH99_Chr5} + 1814975;
$chrtotal{cneoH99_Chr7} = $chrtotal{cneoH99_Chr6} + 1422463;
$chrtotal{cneoH99_Chr8} = $chrtotal{cneoH99_Chr7} + 1399503;
$chrtotal{cneoH99_Chr9} = $chrtotal{cneoH99_Chr8} + 1398693;
$chrtotal{cneoH99_Chr10} = $chrtotal{cneoH99_Chr9} + 1186808;
$chrtotal{cneoH99_Chr11} = $chrtotal{cneoH99_Chr10} + 1059964;
$chrtotal{cneoH99_Chr12} = $chrtotal{cneoH99_Chr11} + 1561994;
$chrtotal{cneoH99_Chr13} = $chrtotal{cneoH99_Chr12} + 774062;
$chrtotal{cneoH99_Chr14} = $chrtotal{cneoH99_Chr13} + 756744;

print OUTPUTALLELES "Chromosome\tPosition\tCum_Position\tCB-100-1\tCB-100-2\tCB-100-3\tCB-100-5\tCB-100-6\tCC-100-1\tCC-100-2\tCC-100-3\tCC-100-4\tCC-100-5\tCC-100-6\tCN-100-1\tCN-100-2\tCN-100-3\tCN-100-4\tCN-100-5\tCN-100-6\t",
"CW-100-1\tCW-100-2\tCW-100-3\tCW-100-4\tCW-100-5\tCW-100-6\tCY-100-1\tCY-100-2\tCY-100-3\tCY-100-4\tCY-100-5\tCY-100-6\tEM3",
"\tP1W10\tP2W10\tP3W10\tP4W10\tP5W10\tP6W10\tSSD719_vitro\tSSD719_vivo\tYSB121\n";

while (<INPUTFILE>)
{
	unless ($_ =~ /#/)
	{
# Gather data
		@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
		$genos[$numstrains+$cols-1] =~ s/\s+$//;
		$chr = $genos[0];
		$chr_name = substr $chr, 11;
		$position = $genos[1];
# Sum up all missing genotypes to determine if locus is low quality
		$missing=0;
		$cum_position = $position + $chrtotal{$chr};
		for ($a=$cols; $a<$cols+$numstrains; $a++)
		{
			@datastrain = split(/:/, $genos[$a]);
			$genotypes[$a]= $datastrain[0];
			$depths[$a] = $datastrain[2];
			$GQ[$a] = $datastrain[3];
			if (($genotypes[$a] =~ /\.\/\./) | ($depths[$a] > $maxdepth) | ($depths[$a] < $mindepth))
			{
				$missing++;
			}
		}
# Need to take care that there could be second variants, indicted by $genos[4] =~ /\,/
		unless ($missing >= $threshold)
		{
			print OUTPUTALLELES "$chr_name\t$position\t$cum_position";
# Next print out the genotypes based on parental and variant columns 3 and 4, respectively
			$refallele = $genos[3];
			if ($genos[4] =~ /\,/)
			{
				@parentals = split(/,/, $genos[4]);
				$parallele1 = $parentals[0];
				$parallele2 = $parentals[1];
				for ($a=$cols; $a<$cols+$numstrains; $a++)
				{
					if (($depths[$a] > $mindepth) && ($GQ[$a] >= $GQmin))
					{
						if ($genotypes[$a] =~ /0\/0/)
						{
							print OUTPUTALLELES "\t$refallele$refallele";
						}
						elsif ($genotypes[$a] =~ /0\/1/)
						{
							print OUTPUTALLELES "\t$refallele$parallele1";
						}
						elsif ($genotypes[$a] =~ /0\/2/)
						{
							print OUTPUTALLELES "\t$refallele$parallele2";
						}	
						elsif ($genotypes[$a] =~ /1\/1/)
						{
							print OUTPUTALLELES "\t$parallele1$parallele1";
						}
						elsif ($genotypes[$a] =~ /1\/2/)
						{
							print OUTPUTALLELES "\t$parallele1$parallele2";
						}
						elsif ($genotypes[$a] =~ /2\/2/)
						{
							print OUTPUTALLELES "\t$parallele2$parallele2";
						}
						else
						{
							print OUTPUTALLELES "\tNA";
						}
					}
					else
					{
						print OUTPUTALLELES "\tNA";
					}
				}
				print OUTPUTALLELES "\n";							
			}
			else
			{
				$parallele = $genos[4];
				for ($a=$cols; $a<$cols+$numstrains; $a++)
				{
					if (($depths[$a] > $mindepth) && ($GQ[$a] >= $GQmin))
					{
						if ($genotypes[$a] =~ /0\/0/)
						{
							print OUTPUTALLELES "\t$refallele$refallele";
						}
						elsif ($genotypes[$a] =~ /0\/1/)
						{
							print OUTPUTALLELES "\t$refallele$parallele";
						}
						elsif ($genotypes[$a] =~ /1\/1/)
						{
							print OUTPUTALLELES "\t$parallele$parallele";
						}
						else
						{
							print OUTPUTALLELES "\tNA";
						}
					}
					else
					{
						print OUTPUTALLELES "\tNA";
					}
				}
				print OUTPUTALLELES "\n";	
			}
		}
	}						
}			
