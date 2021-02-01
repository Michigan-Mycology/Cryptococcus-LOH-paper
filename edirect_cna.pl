#!/usr/bin/perl

# This program takes a list of gene names and searches for entrezids using edirect

# open input and output files
open (INPUTFILE, "gene_names.txt");
open (OUTPUTFILE, "> entrez_ids.txt");

while (<INPUTFILE>)
{
        chomp($_);
        $cmd = 'esearch -db Gene -query ' . $_ . ' | efetch >ids.txt';
        system($cmd);
        open (TEMPFILE, "ids.txt");
        while (<TEMPFILE>)
        {
                print $_;
                if ($_ =~ /ID\:/)
                {
                        $id = substr($_, 4);
                        print OUTPUTFILE "$id";
                }
        }
}


