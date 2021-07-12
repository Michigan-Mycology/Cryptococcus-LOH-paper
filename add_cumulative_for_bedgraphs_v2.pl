#!/usr/bin/perl

# This program takes bedgraph output and adds cumulative length
# It has two genomes and resets after the first genome.

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl add_cumulative_for_bedgraphs.pl inputbedfile outfile\n";
}

# Hash of cumulative lengths
$chrtotal{Crneo_H99_Chr01} = 0;
$chrtotal{Crneo_H99_Chr02} = $chrtotal{Crneo_H99_Chr01} + 2291717;
$chrtotal{Crneo_H99_Chr03} = $chrtotal{Crneo_H99_Chr02} + 1618389;
$chrtotal{Crneo_H99_Chr04} = $chrtotal{Crneo_H99_Chr03} + 1575141;
$chrtotal{Crneo_H99_Chr05} = $chrtotal{Crneo_H99_Chr04} + 1084216;
$chrtotal{Crneo_H99_Chr06} = $chrtotal{Crneo_H99_Chr05} + 1821986;
$chrtotal{Crneo_H99_Chr07} = $chrtotal{Crneo_H99_Chr06} + 1422457;
$chrtotal{Crneo_H99_Chr08} = $chrtotal{Crneo_H99_Chr07} + 1413039;
$chrtotal{Crneo_H99_Chr09} = $chrtotal{Crneo_H99_Chr08} + 1406105;
$chrtotal{Crneo_H99_Chr10} = $chrtotal{Crneo_H99_Chr09} + 1185819;
$chrtotal{Crneo_H99_Chr11} = $chrtotal{Crneo_H99_Chr10} + 1073556;
$chrtotal{Crneo_H99_Chr12} = $chrtotal{Crneo_H99_Chr11} + 1561994;
$chrtotal{Crneo_H99_Chr13} = $chrtotal{Crneo_H99_Chr12} + 774481;
$chrtotal{Crneo_H99_Chr14} = $chrtotal{Crneo_H99_Chr13} + 756939;
$chrtotal{Cryden_JEC21_chrA_1} = 0;
$chrtotal{Cryden_JEC21_chrB_2} = $chrtotal{Cryden_JEC21_chrA_1} + 2300533;
$chrtotal{Cryden_JEC21_chrC_12} = $chrtotal{Cryden_JEC21_chrB_2} + 1632307;
$chrtotal{Cryden_JEC21_chrD_4} = $chrtotal{Cryden_JEC21_chrC_12} + 906719;
$chrtotal{Cryden_JEC21_chrE_5} = $chrtotal{Cryden_JEC21_chrD_4} + 1783081;
$chrtotal{Cryden_JEC21_chrF_6} = $chrtotal{Cryden_JEC21_chrE_5} + 1507550;
$chrtotal{Cryden_JEC21_chrG_7} = $chrtotal{Cryden_JEC21_chrF_6} + 1438950;
$chrtotal{Cryden_JEC21_chrH_9} = $chrtotal{Cryden_JEC21_chrG_7} + 1347793;
$chrtotal{Cryden_JEC21_chrI_10} = $chrtotal{Cryden_JEC21_chrH_9} + 1178688;
$chrtotal{Cryden_JEC21_chrJ_3} = $chrtotal{Cryden_JEC21_chrI_10} + 1085720;
$chrtotal{Cryden_JEC21_chrK_11} = $chrtotal{Cryden_JEC21_chrJ_3} + 2105742;
$chrtotal{Cryden_JEC21_chrL_13} = $chrtotal{Cryden_JEC21_chrK_11} + 1019846;
$chrtotal{Cryden_JEC21_chrM_14} = $chrtotal{Cryden_JEC21_chrL_13} + 787999;
$chrtotal{Cryden_JEC21_chrN_8} = $chrtotal{Cryden_JEC21_chrM_14} + 762694;
# length of JEC21 Chr N(8) = 1194300

print OUTPUTFILE "Chr\tBegin\tEnd\tCum_Position\tDepth\n";

while (<INPUTFILE>)
{
		@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
		$chr = $genos[0];
		$start = $genos[1];
		$end = $genos[2];
		$Depth = $genos[3];
		$cum_mid = $end - (0.5 * ($end - $start)) + $chrtotal{$chr};
		print OUTPUTFILE "$chr\t$start\t$end\t$cum_mid\t$Depth";
}