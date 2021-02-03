#!/usr/bin/perl

use List::Util qw(sum);
use Data::Dumper qw(Dumper);
# use warnings;

# This program takes a file listing breakpoints of LOH and a file listing mutations locations and determines whether these are more spatially clustered than random
# It then generates 1000 random simulated data sets so that the distances can be compared to a random model
# First	argument is a file listing location of mutations across strains
# The second argument is a file listing location of LOH breakpoints across strains

# open input and output files
open (INPUTSNPFILE, $ARGV[0]);
open (INPUTBPFILE, $ARGV[1]);
open (OUTPUTRFILE, "> random_mutation_to_breakpoint_1000X_new.txt");
open (OUTPUTAFILE, "> averages_mutation_to_breakpoint_new.txt");

# Number of random data sets
my $randomizations = 1000;

my $maxbp = 18854554;

# ignore first line
my $header = <INPUTSNPFILE>;
print $header;

# We will sum up lines in SNP file to get number of SNPs needed for randomizations
my $numSNPs = 0;

# This is the hash that stores the location data of mutations
my %snps;

# read in list of snps in columns, number and cumulative location
while (<INPUTSNPFILE>)
{
        chomp($_);
		@line = split(/\t/, $_);
		push @{ $snps{$line[2]} }, $line[1]; 
		$length = scalar( @{ $snps{$line[2]} } );		
		$numSNPs++;
}

# Now read in list of breakpoints and add them to a hash like with the mutations
my $numBPs = 0;

# This is the hash that stores the location data of breakpoints
my %breakpts;

# ignore first line
my $header = <INPUTBPFILE>;
print $header;


while (<INPUTBPFILE>)
{
        chomp($_);
		@line = split(/\t/, $_);
		push @{ $breakpts{$line[2]} }, $line[1]; 
		$length = scalar( @{ $breakpts{$line[2]} } );		
		$numBPs++;
}

# Now go through the SNP keys and calculate a distance from each SNP to each breakpoint
# put that value in array @odist
my @odist = ();

# This number represents the number of SNPs where the distance to a breakpoint is defined.
my $num_measured_SNPs = $numSNPs;

foreach $strain (keys %snps)
{
	foreach $mut_loc ( @{ $snps{$strain} } )
	{
		if ($breakpts{$strain})
		{
# This just sets to the largest value for convenience
			$mindist = $maxbp;
			foreach $ bp_loc ( @{ $breakpts{$strain} } )
			{
				if (abs($mut_loc - $bp_loc) < $mindist)
				{
					$mindist = abs($mut_loc - $bp_loc);
				}
			}
			push (@odist, $mindist);
		}
		else 
		{
			$num_measured_SNPs--;
		}
	}
}

print "number of measured SNPs = $num_measured_SNPs @odist\n";
$mean_real = mean(@odist);
print OUTPUTAFILE "Observed distance from SNP to nearest LOH breakpoint: $mean_real\n";
print "Observed distance from SNP to nearest LOH breakpoint: $mean_real\n";



sub mean {
    return sum(@_)/@_;
}

# Next we will do the randomizations.
# For each of the above, instead of looking at the value of $mut_loc, we will draw a random number between 1 and the genome size.
# Then calculate the distance, average these, and then print to a list of randomizations.

# Store each randomization in an array 
my @randomdists = ();


# Now go through the SNP keys n times where n = # randomizations and calculate a distance from each SNP to each breakpoint

for ($a=0;$a<$randomizations;$a++)
{
# put the temp dist values in array @simdist
	my @simdist = ();
	foreach $strain (keys %snps)
	{
		foreach $mut_loc_2 ( @{ $snps{$strain} } )
		{
# for each mutation choose a new random location in the genome
			my $random_mut = int(rand($maxbp))+1;
			if ($breakpts{$strain})
			{
# This just sets to the largest value for convenience
				$mindist = $maxbp;
				foreach $ bp_loc ( @{ $breakpts{$strain} } )
				{
					if (abs($random_mut - $bp_loc) < $mindist)
					{
						$mindist = abs($random_mut - $bp_loc);
					}
				}
				push (@simdist, $mindist);
			}
		}
	}
	$mean_randomization = mean(@simdist);
	push (@randomdists, $mean_randomization);
	print OUTPUTRFILE $mean_randomization, "\n";
}

$average_randomization = mean(@randomdists);
print OUTPUTAFILE "Mean distance from random SNP to nearest LOH breakpoint: $average_randomization\n";
print "Mean distance from random SNP to nearest LOH breakpoint: $average_randomization\n";