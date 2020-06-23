#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_b $opt_i $opt_e $opt_x $opt_d $opt_a $opt_q $opt_r $opt_f $opt_o $opt_t $opt_k $opt_m);

# Usage
my $usage = "
CombineHetDifsWithFstMatrix.pl - combines the output of HetDifsPerSite.pl and FstMatrix.pl

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

Usage: perl CombineHetDifsWithFstMatrix.pl options
  -b  (path to) outfile baseline name used in HetDifsPerRegion.pl and HetDifsPerSite.pl
  -f  (path to) output of FstPerSite.pl
  -i  (path to) output of FstMatrix.pl
  -e  (path to) output of FreqPerSite.pl
  -d  (path to) chromosome info file
      (tab-delimited table, first column is chrom name, second column is chrom length in bp, third column is centromere start position, fourth column is centromere end position)
  -o  (path to) output file
 optional:
  -x don't inlcude x [default = do include]
  -a don't include autosomes [default = do include]
  -q min frequency to include [default = no min]
  -r max frequency to include [default = no max]
  -t max fst to include [default = no max]
  -k min dif*freq*fstrank to include [default = none]
  -m min number of genotype combinations observed [default = 1000] (this allows filtering of variants with lots of missing data)
";

#############

getopts('b:q:i:e:d:f:o:r:t:m:k:xa');

die $usage unless ($opt_b);
die $usage unless ($opt_f);
die $usage unless ($opt_i);
die $usage unless ($opt_e);
die $usage unless ($opt_d);
die $usage unless ($opt_o);

my ($fst, $matrix, $outgroup, $baselinename, $outfile, $nox, $noauto, $minfreq, $datafile, $maxfreq, $maxfst, $mink, $mincombos);

$fst = $opt_f;

$outfile = $opt_o;

$matrix = $opt_i;

$baselinename = $opt_b;

$datafile = $opt_d if $opt_d;

$outgroup = $opt_e if $opt_e;

if (defined $opt_q) {
  $minfreq = $opt_q;
} else {
  $minfreq = 0;
}

if (defined $opt_m) {
  $mincombos = $opt_m;
} else {
  $mincombos = 1000;
}

if (defined $opt_r) {
  $maxfreq = $opt_r;
} else {
  $maxfreq = 1;
}

if (defined $opt_t) {
  $maxfst = $opt_t;
} else {
  $maxfst = 1;
}

if (defined $opt_k) {
  $mink = $opt_k;
} else {
  $mink = 0;
}

if (defined $opt_x) {
  $nox = 1;
  print "Skipping X\n";
}

if (defined $opt_a) {
  $noauto = 1;
  print "Skipping Autosomes\n";
}

if ((defined $nox)&&(defined $noauto)) {
  print "Can't skip X and autosomes.\n";
  exit;
}

my %chromsizes;

open(DATA, $datafile) || die "can't open $datafile";

while (<DATA>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[1] =~ /\d/) {
    $chromsizes{$data[0]} = $data[1];
  }
}

close (DATA);

my $SNPcount = 0;

my %HoAmatrix;

open(MAT, $matrix) || die "can't open $matrix\n";

while (<MAT>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  my @total = pop @data;
  my $freq = shift @data;
  if ($data[0] =~ /\d/) {
     push @{$HoAmatrix{$freq}}, @data;
  }
}

close (MAT);

unless ( open(OFC, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OFC "SNP\tChrom\tSite\tDif\tFst\tFreqAll\tFstRank\tOutgroupFreq\tAH\tZ\tPause";
close (OFC);

my $extra = 0.01;

foreach my $c (keys %chromsizes) {
  
  my %outgroup;
  
  open(SNPSONE, $outgroup) || die "can't open $outgroup\n";
  
  while (<SNPSONE>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    unless ($data[0] =~ /^$c$/) {
      next;
    }
    if ($data[1] =~ /\d/) {
      my $sitesnp = "$data[1]\t$data[2]";
      $outgroup{$sitesnp} = $data[-1];
    }
  }
  
  close (SNPSONE);
  
  my $file = "$baselinename"."_$c"."_Difs_Per_Site.txt";
  unless(open(SNPS, $file)) {
    print "$file does not exist, skipping.\n";
    next;
  }
  print "Considering $c\n";
  
  my %df;
  
  my %fstfreq;
  
  my $tempsnpcount = 0;

  while (<SNPS>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /\d/) {
      unless (((defined $nox)&&($c =~ /^X$/))||((defined $noauto)&&($c !~ /^X$/))) {
        my $countsum = $data[4] + $data[6] + $data[8];
        if (($countsum >= $mincombos)||(($c =~ /^X$/)&&($countsum >= ($mincombos/2)))) {
          $SNPcount +=1;
          $tempsnpcount +=1;
          my $dif = (int($data[-1]*1000))/1000;
          my $sitesnp = "$data[0]\t$data[1]";
          $df{$sitesnp} = $dif;
        }
      }
    }
  }
  
  close (SNPS);
  
  my @out;
  
  if ($tempsnpcount > 0) {
  
    open(FST, $fst) || die "can't open $fst\n";
    
    while (<FST>) {
      my $line = $_;
      $line =~ s/\r|\n//g;
      my @data = split "\t", $line;
      if ($data[0] =~ /^$c$/) {
        if ($data[-2] =~ /\d/) {
          my $sitesnp = "$data[1]\t$data[2]";
          if (defined $df{$sitesnp}) {
            if (($data[-2] >= $minfreq)&&($data[-2] <= $maxfreq)&&($data[-1] <= $maxfst)) {
              my $freq = int($data[-2]*100);
              if ($freq > 0) {
                my $fsection = "NA";
                for (my $f = 0; $f < 100; $f++) {
                  unless (defined $HoAmatrix{$freq}[$f]) {
                    print "No matrix value for freq $freq position $f!\n$line\n";
                    exit;
                  }
                  if (($data[-1] <= $HoAmatrix{$freq}[$f])&&($data[-1] >= $HoAmatrix{$freq}[$f+1])) {
                    $fsection = $f;
                    last;
                  }
                }
                my $k = $df{$sitesnp}*$data[-2]*($fsection/100);
                if ($k >= $mink) {
                  my $ah = $df{$sitesnp}*$data[-2];
                  unless (defined $outgroup{$sitesnp}) {
                    $outgroup{$sitesnp} = 0;
                  }
                  my $z = ($fsection/100)*($extra/($outgroup{$sitesnp}+$extra));
                  my $pause = $ah*$z;
                  push @out, "$data[2]\t$data[0]\t$data[1]\t$df{$sitesnp}\t$data[-1]\t$data[-2]\t$fsection\t$outgroup{$sitesnp}\t$ah\t$z\t$pause";
                }
              }
            }
          }
        }
      }
    }
    
    close (FST);
  
  }
  
  if (defined $out[0]) {

    my $result = join "\n", @out;
    
    @out = ();
    
    unless ( open(OFC, ">>$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit;
    }
    print OFC "\n$result";
    close (OFC);
  
  }
  
  print "There are $SNPcount total SNPs so far.\n";

}

