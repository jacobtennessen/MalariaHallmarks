#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_f $opt_o $opt_c $opt_a $opt_x);

# Usage
my $usage = "
FstMatrix.pl - Calculates Fst matrix from list of Fsts

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

Usage: perl FstMatrix.pl options
  -f  (path to) FST file produced by FstPerSite.pl
  -o  (path to) outfile for matrix
 optional:
  -c 0-based column for mean frequency [default = -2]
  -x don't inlcude x [default = do include]
  -a don't include autosomes [default = do include]
";

#############

getopts('f:o:c:xa');

die $usage unless ($opt_f);
die $usage unless ($opt_o);

my ($infile, $outfile, $freqcol, $nox, $noauto);

$infile = $opt_f;

$outfile = $opt_o;

if (defined $opt_c) {
  $freqcol = $opt_c;
} else {
  $freqcol = -2;
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

my @out;

my @titlequants;

for (my $q = 0; $q <= 100; $q++) {
  push @titlequants, "FSTq$q";
}

my $titlequants = join "\t", @titlequants;

push @out, "Freq\t$titlequants\tTotal";

my $metatotal = 0;

for (my $f = 1; $f <= 50; $f ++) {

  my @fst;
  
  my $allsnps = 0;
    
  open(IN, $infile) || die "can't open $infile\n";
  
  while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    unless (((defined $nox)&&($data[0] =~ /^X$/))||((defined $noauto)&&($data[0] !~ /^X$/))) {
      if ($data[1] =~ /\d/) {
        my $freq = int($data[$freqcol]*100);
        if ($freq == $f) {
          push @fst, $data[-1];
        }
        $allsnps +=1;
      }
    }
  }
  
  close (IN);
  
  my $total = scalar(@fst);
  
  $metatotal += $total;
  
  print "Addting $total, so seen $metatotal SNPs out of $allsnps at freq $f\n";
  
  if (defined $fst[0]) {
  
    @fst = reverse (sort by_number @fst);

    my @quants;
    
    push @quants, $fst[0];
    
    for (my $q = 1; $q <= 99; $q++) {
      my $pos = sprintf "%.0f", (($q/100)*$total) - 1;
      push @quants, $fst[$pos];
    }
    
    push @quants, $fst[-1];
    
    my $quants = join "\t", @quants;
    
    push @out, "$f\t$quants\t$total";
  
  } else {
    print "Skipping $f\n";
  }

}

my $result = join "\n", @out;

unless ( open(OFC, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OFC $result;
close (OFC);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
