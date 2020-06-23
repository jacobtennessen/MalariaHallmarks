#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_m $opt_a $opt_o $opt_s $opt_c $opt_b $opt_v $opt_d $opt_w $opt_e);

# Usage
my $usage = "
HetDifsPerRegion.pl - get heterozygosity differences from VCF in 1000 Genomes style

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

Usage: perl HetDifsPerRegion.pl options
  -m  directory path for vcfs
  -a  comma-delimited positions of target samples (samples start at 0)
  -v  generic vcf or vcf.gz name with CHROM, e.g. ALL.chrCHROM_GRCh38.genotypes.vcf.gz
                (assumes each chromosome is a separate vcf in metafolder,
                all with the same name expect for the chromosome,
                e.g. ALL.chr1_GRCh38.genotypes.vcf.gz, ALL.chr2_GRCh38.genotypes.vcf.gz, etc.)
  -d  (path to) chromosome info file
      (tab-delimited table, first column is chrom name, second column is chrom length in bp, third column is centromere start position, fourth column is centromere end position)
  -o  (path to) outfile baseline name
 optional:
  -s  size of range to examine [default = 5000000]
  -c  chrom to examine [default = all]
  -b  region to examine [start-end]
  -w  potential distance of a selection signal [default = 1000000]
  -e  distance around centromere to avoid [default = 5000000]
";

#############

getopts('m:a:o:s:c:b:v:d:w:e:');

die $usage unless ($opt_m);
die $usage unless ($opt_a);
die $usage unless ($opt_o);
die $usage unless ($opt_v);
die $usage unless ($opt_d);

my ($metafolder, @pop, $outfilebase, $rangesize, $chromtouse, @targetregion, $vcfname, $datafile, $relevantdist, $centroavoid);

$metafolder = $opt_m if $opt_m;

$outfilebase = $opt_o if $opt_o;

$datafile = $opt_d if $opt_d;

@pop = split ",", $opt_a if $opt_a;

$vcfname = $opt_v;

if (defined $opt_s) {
  $rangesize = $opt_s;
} else {
  $rangesize = 5000000;
}

if (defined $opt_w) {
  $relevantdist = $opt_w;
} else {
  $relevantdist = 1000000;
}

if (defined $opt_e) {
  $centroavoid = $opt_e;
} else {
  $centroavoid = 5000000;
}

if (defined $opt_c) {
  $chromtouse = $opt_c;
}

if (defined $opt_b) {
  @targetregion = split "-", $opt_b;
}

my %chromsizes;

my %centromeres;

open(DATA, $datafile) || die "can't open $datafile";

while (<DATA>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[1] =~ /\d/) {
    $chromsizes{$data[0]} = $data[1];
    my $centromere = "$data[2]-$data[3]";
    $centromeres{$data[0]} = $centromere;
  }
}

close (DATA);

foreach my $chrom (keys %chromsizes) {
  
  if ((defined $chromtouse)&&($chromtouse !~ /^$chrom$/)) {
    next;
  }
  
  my $chromvcf = $vcfname;
  
  $chromvcf =~ s/CHROM/$chrom/g;
  
  my $genotypes = "$metafolder$chromvcf";
  
  my @fullcentromere = split "-", $centromeres{$chrom};
  
  push @fullcentromere, $chromsizes{$chrom};
  
  unshift @fullcentromere, 0;
  
  my @regions;
  
  for (my $f = 0; $f < (scalar(@fullcentromere)); $f+=2) {
      
    my $regstart = $fullcentromere[$f];
    
    unless ($regstart == 0) {
      $regstart += $centroavoid;
    }
    
    my $regend = $fullcentromere[$f+1];
    
    unless ($regend == $chromsizes{$chrom}) {
      $regend -= $centroavoid;
    }
    
    for (my $r = $regstart; $r < ($regend - $rangesize); $r += ($rangesize - 2*$relevantdist)) {
      if (defined $targetregion[0]) {
        if ((($r + $rangesize - $relevantdist) < $targetregion[0])||(($r + $relevantdist) > $targetregion[1])) {
          next;
        }
      }
      push @regions, $r;
    }
    
    my $closestbefore = ($regend - $rangesize);
    
    unless ($regions[-1] >= ($closestbefore - $relevantdist)) {
      push @regions, $closestbefore;
    }
  
  }
  
  foreach my $r (@regions) {
  
    my @allsamples;
    
    if ($genotypes =~ /.gz$/) {
      open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
    }
    else {
      open(IN, $genotypes) || die "can't open $genotypes";
    }
    
    my %hets;
    
    while (<IN>) {
      my $line = $_;
      $line =~ s/\r|\n//g;
      my @data = split "\t", $line;
      if ($data[0] =~ /^#CHROM/) {
        for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
            push @allsamples, $fnx;
        }
        next;
      } elsif ($line =~ /^#/) {
          next;
      }
      if (($data[1] < $r)||($data[1] > ($r + $rangesize))) {
        next;
      }
      my @inds = @data[@allsamples[@pop]];
      my $flag = 0;
      my $seenhemi = 0;
      if ($chrom =~ /\d/) {
        foreach my $i1 (@inds) {
          unless ($i1 =~ /\d\|\d/) {
            $flag = 1;
            last;
          }
        }
      } else {
        foreach my $i1 (@inds) {
          unless (($i1 =~ /^\d$/)||($i1 =~/\d\|\d/)) {
            $flag = 1;
            last;
          }
          if ($i1 =~ /^\d$/) {
            $seenhemi = 1;
          }
        }
      }
      if (($chrom !~ /\d/)&&($seenhemi == 0)) {
        $flag = 1;
      }
      if ($flag == 0) {
        my @seeninds;
        my $ind1count = 0;
        foreach my $i1 (@inds) {
          my @i1data = split /\|/, $i1;
          my $ind2count = 0;
          if (defined $i1data[1]) {
            unless ($i1data[0] == $i1data[1]) {
              my $combo = "$pop[$ind1count]"."_1\t$pop[$ind1count]"."_2";
              if (defined $hets{$combo}) {
                $hets{$combo} += 1;
              } else {
                $hets{$combo} = 1;
              }
            }
          }
          foreach my $i2 (@seeninds) {
            my @i2data = split /\|/, $i2;
            unless ($i2data[0] == $i1data[0]) {
              my $combo = "$pop[$ind2count]"."_1\t$pop[$ind1count]"."_1";
              if (defined $hets{$combo}) {
                $hets{$combo} += 1;
              } else {
                $hets{$combo} = 1;
              }
            }
            if (defined $i1data[1]) {
              unless ($i2data[0] == $i1data[1]) {
                my $combo = "$pop[$ind2count]"."_1\t$pop[$ind1count]"."_2";
                if (defined $hets{$combo}) {
                  $hets{$combo} += 1;
                } else {
                  $hets{$combo} = 1;
                }
              }
              if (defined $i2data[1]) {
                unless ($i2data[1] == $i1data[1]) {
                  my $combo = "$pop[$ind2count]"."_2\t$pop[$ind1count]"."_2";
                  if (defined $hets{$combo}) {
                    $hets{$combo} += 1;
                  } else {
                    $hets{$combo} = 1;
                  }
                }
              }
            }
            if (defined $i2data[1]) {
              unless ($i2data[1] == $i1data[0]) {
                my $combo = "$pop[$ind2count]"."_2\t$pop[$ind1count]"."_1";
                if (defined $hets{$combo}) {
                  $hets{$combo} += 1;
                } else {
                  $hets{$combo} = 1;
                }
              }
            }
            $ind2count +=1;
          }
          push @seeninds, $i1;
          $ind1count +=1;
        }
      }
    }
    
    close (IN);
    
    my @out;
    
    for (my $p1 = 0; $p1 < (scalar(@pop)); $p1++) {
      for (my $p2 = $p1; $p2 < (scalar(@pop)); $p2++) {
        if ($p1 == $p2) {
          my $comboself = "$pop[$p1]"."_1\t$pop[$p2]"."_2";
          unless (defined $hets{$comboself}) {
            $hets{$comboself} = 0;
          }
          push @out, "$comboself\t$hets{$comboself}";
        } else {
          my $combo11 = "$pop[$p1]"."_1\t$pop[$p2]"."_1";
          unless (defined $hets{$combo11}) {
            $hets{$combo11} = 0;
          }
          my $combo12 = "$pop[$p1]"."_1\t$pop[$p2]"."_2";
          unless (defined $hets{$combo12}) {
            $hets{$combo12} = 0;
          }
          my $combo21 = "$pop[$p1]"."_2\t$pop[$p2]"."_1";
          unless (defined $hets{$combo21}) {
            $hets{$combo21} = 0;
          }
          my $combo22 = "$pop[$p1]"."_2\t$pop[$p2]"."_2";
          unless (defined $hets{$combo22}) {
            $hets{$combo22} = 0;
          }
          push @out, "$combo11\t$hets{$combo11}";
          push @out, "$combo12\t$hets{$combo12}";
          push @out, "$combo21\t$hets{$combo21}";
          push @out, "$combo22\t$hets{$combo22}";
        }
      }
    }
    
    my $result = join "\n", @out;
    
    @out = ();
    
    my $end = $r + $rangesize;
    
    my $outfile = "$outfilebase"."_$chrom"."_$r-$end.txt";
    
    unless ( open(SUMM, ">$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit
    }
    print SUMM "Ind1\tInd2\tHet\n$result";
    close (SUMM);
  
  }
}
