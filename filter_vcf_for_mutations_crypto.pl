#!/usr/bin/perl

# This program takes a vcf file and filters it on quality to identify novel mutations
# Qualities
# 1. All ancestors need to be called and 99 GQ and 0/0 and depth at least 18X
# 2. Any mutation needs to be GQ 99, 0/1, or 1/1, and depth at least 18X and not greater than 350X
# 3. No missing data for any strain. i.e., No ./. 

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

if (@ARGV < 1) {
    die "Improper number of arguments. Usage: perl filter_vcf_for_mutations_crypto.pl inputvcffile outputvcffile\n";
}

$numstrains=39;

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# mindepth is the minimum number of reads for a strain to be retained
$mindepth = 16;

# maxdepthpool is the maximum number of reads or it should be filtered out
$maxdepth = 350;

# Threshold is the max number of strains for which missing data is possible for a locus to be included
$threshold = 1;

# Set this to have a minimum quality of a genotype to be used. The GQ score of GATK
$GQmin = 99;

# End of evol strains part 1
$endstrain1 = 37;

# End of evol strains part 2
$endstrain2 = 44;

# Column for YSB121
$Aparentcol = 47;

# Column for EM3
$Dparentcol = 38;

# Column for SSD719-vitro
$dipvitroparentcol = 45;

# Column for SSD719-vivo
$dipvivoparentcol =46;

while (<INPUTFILE>)
{
	if ($_ =~ /#/)
	{
		print OUTPUTFILE $_;
	}
	else
	{
# Splits up the line of data by tabs
		@genos = split(/\t/, $_);
#	print "@genos\n";
		$chr = $genos[0];
		$chr_name = substr $chr, 11;
		$position = $genos[1];
		print "\n$chr\t$position";
# This following line removes the whitespace at the end of the line.
		$genos[$numstrains+$cols-1] =~ s/\s+$//;
# Sum up all missing genotypes to determine if locus is low quality
		$missing=0;
# Indicates one or more parents failed tests 
		$badparent = "false";
# Flag for found a mutant
		$foundmutant = "false";
		for ($a=$cols; $a<$cols+$numstrains; $a++)
		{
			@datastrain = split(/:/, $genos[$a]);
			$genotypes[$a]= $datastrain[0];
			$depths[$a] = $datastrain[2];
			$GQ[$a] = $datastrain[3];
			if ($genotypes[$a] =~ /\.\/\./) 
			{
				$missing++;
			}
		}
# Test parent genotype and qual
		if (($genotypes[$Dparentcol] =~ /1/) | ($genotypes[$Aparentcol] =~ /1/) | ($genotypes[$dipvitroparentcol] =~ /1/) | ($genotypes[$dipvivoparentcol] =~ /1/))
		{
			$badparent = "true";
		}
		if (($depths[$Dparentcol] < $mindepth) | ($depths[$Aparentcol] < $mindepth) | ($depths[$dipvitroparentcol] < $mindepth) | ($depths[$dipvivoparentcol] < $mindepth))
		{
			$badparent = "true";

		}
		if (($GQ[$Dparentcol] < $GQmin) | ($GQ[$Aparentcol] < $GQmin) | ($GQ[$dipvitroparentcol] < $GQmin) | ($GQ[$dipvivoparentcol] < $GQmin))
		{
			$badparent = "true";
		}			
# Now check that there is one strain with mutation which is high quality
		unless (($missing >= $threshold) | ($badparent eq "true"))
		{ 
			for ($a=$cols; $a<$cols+$endstrain1; $a++)
			{
				if (($genotypes[$a] =~ /1/) && ($GQ[$a] >= $GQmin) && ($depths[$a] >= $mindepth))
				{
					$foundmutant = "true";
				}
			}	
			for ($a=$Dparentcol+1; $a<$cols+$endstrain2; $a++)
			{
				if (($genotypes[$a] =~ /1/) && ($GQ[$a] >= $GQmin) && ($depths[$a] >= $mindepth))
				{
					$foundmutant = "true";
				}
			}	

			if ($foundmutant eq "true")
			{
				print OUTPUTFILE $_;
			}
		}	
	}						
}			
