#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_o $opt_a $opt_b $opt_c $opt_r $opt_v $opt_x $opt_m);

# Usage
my $usage = "
GetTrident.pl - Finds windows with high 3-way Fst

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

Usage: perl GetTrident.pl options
  -a  (path to) first Fst region file
  -b  (path to) second Fst region file
  -c  (path to) third Fst region file
  -v  (path to) list of viable regions
  -o  (path to) outfile
  -r  region size
 optional:
  -x  include X chromosome [default = no]
  -m  mimimum fst to include [default = -1]
";

#############

getopts('o:a:b:c:r:v:m:x');

die $usage unless ($opt_o);
die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_c);
die $usage unless ($opt_r);

my ($outfile, $fst1, $fst2, $fst3, $viable, $regionsize, $includex, $minfst);

$outfile = $opt_o if $opt_o;

$fst1 = $opt_a if $opt_a;

$fst2 = $opt_b if $opt_b;

$fst3 = $opt_c if $opt_c;

$viable = $opt_v if $opt_v;

$regionsize = $opt_r if $opt_r;

if (defined $opt_x) {
  $includex = 1;
};

if (defined $opt_m) {
  $minfst = $opt_m;
} else {
  $minfst = -1;
}

my %HoHviable;

my $totalwindows = 0;

open(VIAB, $viable) || die "can't open $viable\n";

while (<VIAB>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[1] =~ /\d/) {
    if ((defined $includex)||($data[0] =~ /\d/)) {
      my $region;
      my $end = $data[1] + $regionsize;
      $region = "$data[1]-$end";
      $HoHviable{$data[0]}{$region} = 1;
      $totalwindows +=1;
    }
  }
}

close (VIAB);

print "Found $totalwindows viable regions.\n";

my %fst1;

my @fst1;

my %pos;

open(SNPSONE, $fst1) || die "can't open $fst1\n";

while (<SNPSONE>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  my $region = "$data[2]-$data[3]";
  if (($data[4] =~ /\d/)&&($data[4] >= $minfst)&&(defined $HoHviable{$data[1]}{$region})) {
    $fst1{$data[0]} = $data[4];
    push @fst1, $data[4];
    $pos{$data[0]} = "$data[1]\t$data[2]\t$data[3]";
  }
}

close (SNPSONE);

@fst1 = reverse (sort by_number @fst1);

my $rank = 0;

my %ranks1;

foreach my $fst1 (@fst1) {
  $rank +=1;
  $ranks1{$fst1} = $rank;
}

my %windowscore1;

foreach my $window1 (keys %fst1) {
  $windowscore1{$window1} = $ranks1{$fst1{$window1}}/$totalwindows;
}

my %fst2;

my @fst2;

open(SNPSTWO, $fst2) || die "can't open $fst2\n";

while (<SNPSTWO>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  my $region = "$data[2]-$data[3]";
  if (($data[4] =~ /\d/)&&($data[4] >= $minfst)&&(defined $HoHviable{$data[1]}{$region})) {
    $fst2{$data[0]} = $data[4];
    push @fst2, $data[4];
  }
}

close (SNPSTWO);

@fst2 = reverse (sort by_number @fst2);

$rank = 0;

my %ranks2;

foreach my $fst2 (@fst2) {
  $rank +=1;
  $ranks2{$fst2} = $rank;
}

my %windowscore2;

foreach my $window2 (keys %fst2) {
  $windowscore2{$window2} = $ranks2{$fst2{$window2}}/$totalwindows;
}

my %fst3;

my @fst3;

open(SNPSTHREE, $fst3) || die "can't open $fst3\n";

while (<SNPSTHREE>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  my $region = "$data[2]-$data[3]";
  if (($data[4] =~ /\d/)&&($data[4] >= $minfst)&&(defined $HoHviable{$data[1]}{$region})) {
    $fst3{$data[0]} = $data[4];
    push @fst3, $data[4];
  }
}

close (SNPSTHREE);

@fst3 = reverse (sort by_number @fst3);

$rank = 0;

my %ranks3;

foreach my $fst3 (@fst3) {
  $rank +=1;
  $ranks3{$fst3} = $rank;
}

my %windowscore3;

foreach my $window3 (keys %fst3) {
  $windowscore3{$window3} = $ranks3{$fst3{$window3}}/$totalwindows;
}

my @out;

foreach my $g (keys %windowscore1) {
  if ((defined $windowscore2{$g})&&(defined $windowscore3{$g})) {
    my $lowest = $windowscore1{$g};
    if ($windowscore2{$g} < $lowest) {
      $lowest = $windowscore2{$g};
    }
    if ($windowscore3{$g} < $lowest) {
      $lowest = $windowscore3{$g};
    }
    my $highest = $windowscore1{$g};
    if ($windowscore2{$g} > $highest) {
      $highest = $windowscore2{$g};
    }
    if ($windowscore3{$g} > $highest) {
      $highest = $windowscore3{$g};
    }
    my $trident = (($lowest+$highest)/2)*(($lowest+$highest)/2);
    my $corrected = sprintf "%.7f", $trident*$totalwindows;
    push @out, "$g\t$pos{$g}\t$fst1{$g}\t$fst2{$g}\t$fst3{$g}\t$ranks1{$fst1{$g}}\t$ranks2{$fst2{$g}}\t$ranks3{$fst3{$g}}\t$trident\t$corrected";
  }
}

my $result = join "\n", @out;

@out = ();

unless ( open(SUMM, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print SUMM "Window\tChrom\tStart\tEnd\tFst1\tFst2\tFst3\tRank1\tRank2\tRank3\tTrident\tCorrectedTrident\n$result";
close (SUMM);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
