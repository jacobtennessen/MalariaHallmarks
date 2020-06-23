#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_m $opt_v $opt_o $opt_p $opt_q $opt_r $opt_l);

# Usage
my $usage = "
FreqPerSite.pl - Calculates minor allele frequency from original vcfs

Copyright (C) 2020 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl FreqPerSite.pl options
  -m  metafolder
  -v  generic vcf or vcf.gz name with CHROM, e.g. ALL.chrCHROM_GRCh38.genotypes.vcf.gz
                (assumes each chromosome is a separate vcf in metafolder,
                all with the same name expect for the chromosome,
                e.g. ALL.chr1_GRCh38.genotypes.vcf.gz, ALL.chr2_GRCh38.genotypes.vcf.gz, etc.)
  -p  comma-delimited lists of positions (samples start at 0) (e.g. 0,3,4,7,8,9)
  -o  (path to) outfile
 optional:
  -q  minimum MAF for export to candidate list [default > 0]
  -r  maximum MAF for export to candidate list [default = 0.5]
  -l  comma-delimited list of chromosomes to examine (default = X,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
";

#############

getopts('m:v:o:p:q:r:l:');

die $usage unless ($opt_m);
die $usage unless ($opt_o);
die $usage unless ($opt_p);
die $usage unless ($opt_v);

my ($metafolder, $outfile, @pop, $minfreq, $maxfreq, $vcfname, @chroms);

$metafolder = $opt_m;

$outfile = $opt_o;

@pop = split ",", $opt_p;

$vcfname = $opt_v;

if (defined $opt_l) {
  @chroms = split ",", $opt_l;
} else {
  @chroms = ("X",1..22);
}

if (defined $opt_q) {
  $minfreq = $opt_q;
} else {
  $minfreq = 0.00000001;
}

if (defined $opt_r) {
  $maxfreq = $opt_r;
} else {
  $maxfreq = 0.5;
}

my @allsamples;

my $seensnps = 0;

unless ( open(OFC, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OFC "Chrom\tSite\tSNP\tRef\tAlt\tSize\tFreq";
close (OFC);

foreach my $c (@chroms) {
  
  my @outliers;

  my $chromvcf = $vcfname;
  
  $chromvcf =~ s/CHROM/$c/g;
  
  my $genotypes = "$metafolder$chromvcf";
  
  if ($genotypes =~ /.gz$/) {
    open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
  }
  else {
    open(IN, $genotypes) || die "can't open $genotypes";
  }
  
  while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /^#CHROM/) {
      unless (defined $allsamples[0]) {
        for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
            push @allsamples, $fnx;
        }
      }
      next;
    } elsif ($line =~ /^#/) {
        next;
    }
    my @inds1 = @data[@allsamples[@pop]];
    my $genoflag = 0;
    my $seenref = 0;
    my $seenalt = 0;
    foreach my $i1 (@inds1) {
      if ($i1 =~ /0\|0/) {
        $seenref += 2;
      } elsif (($i1 =~ /0\|1/)||($i1 =~ /1\|0/)||($i1 =~ /2\|0/)||($i1 =~ /0\|2/)) {
        $seenref += 1;
        $seenalt += 1;
      } elsif (($i1 =~ /1\|1/)||($i1 =~ /1\|2/)||($i1 =~ /2\|1/)||($i1 =~ /2\|2/)) {
        $seenalt += 2;
      } elsif ($i1 =~ /^0$/) {
        $seenref +=1;
      } elsif (($i1 =~ /^1$/)||($i1 =~ /^2$/)) {
        $seenalt +=1;
      } else {
        $genoflag += 1;
      }
    }
    if (($genoflag == 0)&&($seenref > 0)&&($seenalt > 0)) {
      my $size = $seenalt+$seenref;
      my $freq = sprintf "%.4f", $seenalt/$size;
      if ($freq > 0.5) {
        $freq = 1 - $freq;
      }
      if (($freq >= $minfreq)&&($freq <= $maxfreq)) {
        push @outliers, "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$size\t$freq";
        $seensnps +=1;
      }
    }
  }
  
  close (IN);
  
  my $outliers = join "\n", @outliers;
  
  @outliers = ();
  
  unless ( open(OFC, ">>$outfile") ) {
      print "Cannot open file \"$outfile\" to write to!!\n\n";
      exit;
  }
  print OFC "\n$outliers";
  close (OFC);
  
  print "Chrom $c complete; $seensnps SNPs seen.\n";

}
