#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_m $opt_o $opt_a $opt_b $opt_c $opt_f $opt_r $opt_s $opt_v $opt_l);

# Usage
my $usage = "
FindViableFstRegionsForTrident.pl - finds regions that could potentially show high Fst among three populations

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

Usage: perl FindViableFstRegionsForTrident.pl options
  -m  directory path for vcfs
  -a  comma-delimited positions, first pop (samples start at 0)
  -b  comma-delimited positions, second pop (samples start at 0)
  -c  comma-delimited positions, third pop (samples start at 0)
  -v  generic vcf or vcf.gz name with CHROM, e.g. ALL.chrCHROM_GRCh38.genotypes.vcf.gz
                (assumes each chromosome is a separate vcf in metafolder,
                all with the same name expect for the chromosome,
                e.g. ALL.chr1_GRCh38.genotypes.vcf.gz, ALL.chr2_GRCh38.genotypes.vcf.gz, etc.)
  -o  (path to) outfile for list of regions
 optional:
  -f  minimum freq to be viable [default = 5%]
  -r  region size [default = 5000]
  -s  stagger offset [default = 0]
  -l  comma-delimited list of chromosomes to examine (default = X,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
";

#############

getopts('m:o:a:b:c:f:r:s:v:l:');

die $usage unless ($opt_m);
die $usage unless ($opt_o);
die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_c);
die $usage unless ($opt_v);

my ($metafolder, $outfile, @pop1, @pop2, @pop3, $minfreq, $regionsize, $offset, $vcfname, @chroms);

$metafolder = $opt_m;

$outfile = $opt_o;

@pop1 = split ",", $opt_a;

@pop2 = split ",", $opt_b;

@pop3 = split ",", $opt_c;

$vcfname = $opt_v;

if (defined $opt_f) {
  $minfreq = $opt_f;
} else {
  $minfreq = 0.05;
}

if (defined $opt_r) {
  $regionsize = $opt_r;
} else {
  $regionsize = 5000;
}

if (defined $opt_s) {
  $offset = $opt_s;
} else {
  $offset = 0;
}

if (defined $opt_l) {
  @chroms = split ",", $opt_l;
} else {
  @chroms = ("X",1..22);
}

my $namelist1;

my $namelist2;

my $namelist3;

my @allsamples;

unless ( open(OFC, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OFC "Chrom\tStart\tSNPs";
close (OFC);

foreach my $c (@chroms) {
  
  my %snps;
  
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
      unless (defined $namelist1) {
        for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
            push @allsamples, $fnx;
        }
        my @names = @data[@allsamples];
        $namelist1 = join ",", @names[@pop1];
        $namelist2 = join ",", @names[@pop2];
        $namelist3 = join ",", @names[@pop3];
      }
      next;
    } elsif ($line =~ /^#/) {
        next;
    }
    my @inds1 = @data[@allsamples[@pop1]];
    my @inds2 = @data[@allsamples[@pop2]];
    my @inds3 = @data[@allsamples[@pop3]];
    my $genoflag = 0;
    my $seenref1 = 0;
    my $seenalt1 = 0;
    my $seenref2 = 0;
    my $seenalt2 = 0;
    my $seenref3 = 0;
    my $seenalt3 = 0;
    foreach my $i1 (@inds1) {
      if ($i1 =~ /0\|0/) {
        $seenref1 += 2;
      } elsif (($i1 =~ /0\|1/)||($i1 =~ /1\|0/)||($i1 =~ /2\|0/)||($i1 =~ /0\|2/)) {
        $seenref1 += 1;
        $seenalt1 += 1;
      } elsif (($i1 =~ /1\|1/)||($i1 =~ /1\|2/)||($i1 =~ /2\|1/)||($i1 =~ /2\|2/)) {
        $seenalt1 += 2;
      } elsif ($i1 =~ /^0$/) {
        $seenref1 +=1;
      } elsif (($i1 =~ /^1$/)||($i1 =~ /^2$/)) {
        $seenalt1 +=1;
      } else {
        $genoflag += 1;
      }
    }
    foreach my $i2 (@inds2) {
      if ($i2 =~ /0\|0/) {
        $seenref2 += 2;
      } elsif (($i2 =~ /0\|1/)||($i2 =~ /1\|0/)||($i2 =~ /2\|0/)||($i2 =~ /0\|2/)) {
        $seenref2 += 1;
        $seenalt2 += 1;
      } elsif (($i2 =~ /1\|1/)||($i2 =~ /1\|2/)||($i2 =~ /2\|1/)||($i2 =~ /2\|2/)) {
        $seenalt2 += 2;
      } elsif ($i2 =~ /^0$/) {
        $seenref2 +=1;
      } elsif (($i2 =~ /^1$/)||($i2 =~ /^2$/)) {
        $seenalt2 +=1;
      } else {
        $genoflag += 1;
      }
    }
    foreach my $i3 (@inds3) {
      if ($i3 =~ /0\|0/) {
        $seenref3 += 2;
      } elsif (($i3 =~ /0\|1/)||($i3 =~ /1\|0/)||($i3 =~ /2\|0/)||($i3 =~ /0\|2/)) {
        $seenref3 += 1;
        $seenalt3 += 1;
      } elsif (($i3 =~ /1\|1/)||($i3 =~ /1\|2/)||($i3 =~ /2\|1/)||($i3 =~ /2\|2/)) {
        $seenalt3 += 2;
      } elsif ($i3 =~ /^0$/) {
        $seenref3 +=1;
      } elsif (($i3 =~ /^1$/)||($i3 =~ /^2$/)) {
        $seenalt3 +=1;
      } else {
        $genoflag += 1;
      }
    }
    if (($genoflag == 0)&&(($seenref1+$seenref2+$seenref3) > 0)&&(($seenalt1+$seenalt2+$seenalt3) > 0)) {
      my $size1 = $seenalt1+$seenref1;
      my $freq1 = $seenalt1/$size1;
      my $size2 = $seenalt2+$seenref2;
      my $freq2 = $seenalt2/$size2;
      my $size3 = $seenalt3+$seenref3;
      my $freq3 = $seenalt3/$size3;
      my $meanfreq = ($freq1+$freq2+$freq3)/3;
      if (($meanfreq >= $minfreq)&&($meanfreq <= (1-$minfreq))) {
        my $region = $regionsize*int(($data[1]-$offset)/$regionsize)+$offset;
        if (defined $snps{$region}) {
          $snps{$region} += 1;
        } else {
          $snps{$region} = 1;
        }
      }
    }
  }
  
  close (IN);
  
  my @goodregions;
  
  foreach my $r (keys %snps) {
    if ($snps{$r} >= 2) {
      push @goodregions, "$c\t$r\t$snps{$r}";
    }
  }
  
  my $goodregions = join "\n", @goodregions;
  
  @goodregions = ();
  
  unless ( open(OFC, ">>$outfile") ) {
      print "Cannot open file \"$outfile\" to write to!!\n\n";
      exit;
  }
  print OFC "\n$goodregions";
  close (OFC);

}

