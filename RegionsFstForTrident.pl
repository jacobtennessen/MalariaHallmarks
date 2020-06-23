#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_o $opt_a $opt_x $opt_m $opt_j $opt_k $opt_l $opt_r $opt_s);

# Usage
my $usage = "
RegionsFstForTrident.pl - Finds highest Fst SNP for regions

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

Usage: perl RegionsFstForTrident.pl options
  -a  (path to) Fst file (output of FstForTrident.pl)
  -o  (path to) outfile
 optional:
  -m  minimum Fst [default = 0.1]
  -j  0-delimted chromsite [default = 0]
  -k  0-delimted snp site [default = 1]
  -l  0-delimted fst site [default = 10]
  -r  region size [default = 5000]
  -s  stagger offset [default = 0]
  -x  include X chrom
";

#############

getopts('o:a:m:j:k:l:r:s:x');

die $usage unless ($opt_o);
die $usage unless ($opt_a);

my ($outfile, $fst, $usex, $minfst, $chromsite, $snpsite, $fstsite, $regionsize, $offset);

$outfile = $opt_o;

$fst = $opt_a if $opt_a;

if (defined $opt_x) {
  $usex = 1;
}

if (defined $opt_m) {
  $minfst = $opt_m;
} else {
  $minfst = 0.1;
}

if (defined $opt_j) {
  $chromsite = $opt_j;
} else {
  $chromsite = 0;
}

if (defined $opt_k) {
  $snpsite = $opt_k;
} else {
  $snpsite = 1;
}

if (defined $opt_l) {
  $fstsite = $opt_l;
} else {
  $fstsite = 10;
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

my %HoHregions;

my $origcount = 0;

open(SNPSONE, $fst) || die "can't open $fst\n";

while (<SNPSONE>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  unless ((defined $usex)||($data[0] =~ /\d/)) {
    next;
  }
  if (($data[$snpsite] =~ /\d/)&&($data[$fstsite] >= $minfst)) {
    my $region = $regionsize*int(($data[$snpsite]-$offset)/$regionsize)+$offset;
    unless ((defined $HoHregions{$data[$chromsite]}{$region})&&($HoHregions{$data[$chromsite]}{$region} >= $data[$fstsite])) {
      $HoHregions{$data[$chromsite]}{$region} = $data[$fstsite];
    }
    $origcount +=1;
  }
}

close (SNPSONE);

print "Found $origcount SNPs in $fst\n";

my @out;

foreach my $chrom (keys %HoHregions) {
  foreach my $region (keys %{ $HoHregions{$chrom} }) {
    my $name = "$chrom"."_$region";
    my $end = $region + $regionsize;
    push @out, "$name\t$chrom\t$region\t$end\t$HoHregions{$chrom}{$region}\t$region\t$regionsize"
  }
}

my $result = join "\n", @out;

@out = ();

unless ( open(SUMM, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print SUMM "Window\tChrom\tStart\tEnd\tFst\tSite\tLength\n$result";
close (SUMM);


