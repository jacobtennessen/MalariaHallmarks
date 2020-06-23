#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_v $opt_a $opt_w $opt_o);

# Usage
my $usage = "
GetDango.pl - calculates sum of linkage disequilibrium scores within a given window for each variant in a vcf

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

Usage: perl GetDango.pl options
  -v  (path to) vcf (can be .gz)
  -a  comma-delimited positions of target samples (samples start at 0)
  -o  (path to) outfile
   optional:
  -w  window size to examine [default = 500]
";

#############

getopts('v:a:w:o:');

die $usage unless ($opt_v);
die $usage unless ($opt_a);
die $usage unless ($opt_o);

my ($genotypes, @pop, $window, $outfile);

$genotypes = $opt_v;

$outfile = $opt_o;

@pop = split ",", $opt_a;

if (defined $opt_w) {
    $window = $opt_w;
} else {
    $window = 500;
}

my @gametelists;

my @sites;

my @rs;

my @deltasum;

my @allsamples;

my @buddies;

my $retainedcount = 0;

my $chrom;

my $gametecountauto = 2*(scalar(@pop));

my $gametecountx;

my @out;

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
    for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
        push @allsamples, $fnx;
    }
    my @names = @data[@allsamples];
    my $namelist = join ",", @names[@pop];
    unless ( open(SUMM, ">$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit;
    }
    print SUMM "#$namelist\n";
    print SUMM "Chrom\tSite\tSNP\tDango\tStotal";
    close (SUMM);
    next;
  } elsif ($line =~ /^#/) {
      next;
  }
  unless (defined $chrom) {
    $chrom = $data[0];
  }
  my $discard = 0;
  for (my $s = 0; $s < (scalar(@sites)); $s++) {
    if ($sites[$s] < ($data[1] - $window)) {
      $discard +=1;
    }
  }
  if ($discard > 0) {
    for (my $d = 1; $d <= $discard; $d++) {
      my $tempsize = scalar(@gametelists);
      my $deltasum = sprintf "%.4f", $deltasum[0];
      push @out, "$chrom\t$sites[0]\t$rs[0]\t$deltasum\t$buddies[0]";
      my $discardlist = shift @gametelists;
      my $discardsite = shift @sites;
      my $discardrs = shift @rs;
      my $discarddeltasum = shift @deltasum;
      my $discardbuddies = shift @buddies;
      $retainedcount -=1;
    }
  }
  my @inds = @data[@allsamples[@pop]];
  my %allelecounts;
  my @templist;
  my $flag = 0;
  my $seenhemi = 0;
  if ($chrom =~ /\d/) {
    foreach my $i (@inds) {
      unless ($i =~ /\d\|\d/) {
        $flag = 1;
        last;
      }
      my @alleles = split /\|/, $i;
      foreach my $a (@alleles) {
        push @templist, $a;
        if (defined $allelecounts{$a}) {
          $allelecounts{$a} += 1;
        } else {
          $allelecounts{$a} = 1;
        }
      }
    }
  } else {
    my $gcx = 0;
    foreach my $i (@inds) {
      unless (($i =~ /^\d$/)||($i =~/\d\|\d/)) {
        $flag = 1;
        last;
      }
      if ($i =~ /^\d$/) {
        $seenhemi = 1;
        $gcx +=1;
      } else {
        $gcx +=2;
      }
      my @alleles = split /\|/, $i;
      foreach my $a (@alleles) {
        push @templist, $a;
        if (defined $allelecounts{$a}) {
          $allelecounts{$a} += 1;
        } else {
          $allelecounts{$a} = 1;
        }
      }
    }
    if ($seenhemi == 0) {
      $flag = 1;
    } else {
      unless (defined $gametecountx) {
        $gametecountx = $gcx;
      }
    }
  }
  if (($flag == 0)&&(scalar(keys %allelecounts) > 1)) {
    my $major;
    my $majorcount = 0;
    foreach my $a (keys %allelecounts) {
      if ($allelecounts{$a} > $majorcount) {
        $majorcount = $allelecounts{$a};
        $major = $a;
      }
    }
    my @finallist;
    foreach my $t (@templist) {
      if ($t == $major) {
        push @finallist, 0;
      } else {
        push @finallist, 1;
      }
    }
    my $finallist = join ",", @finallist;
    push @gametelists, $finallist;
    push @sites, $data[1];
    push @rs, $data[2];
    push @deltasum, 0;
    push @buddies, $retainedcount;
    $retainedcount += 1;
    foreach my $b (@buddies) {
      $b +=1;
    }
    my $gametecount = $gametecountauto;
    if ($chrom !~ /\d/) {
       $gametecount = $gametecountx;
    }
    for (my $si = 1; $si < $retainedcount; $si++) {
      my $ci = 0;
      my $cj = 0;
      my $cij = 0;
      my @gametesi = split ",", $gametelists[$si-1];
      for (my $n = 0; $n < $gametecount; $n++) {
        $ci += $gametesi[$n];
        $cj += $finallist[$n];
        if (($gametesi[$n] == 1)&&($finallist[$n] == 1)) {
          $cij +=1;
        }
      }
      my $pi = $ci/$gametecount;
      my $pj = $cj/$gametecount;
      my $pij = $cij/$gametecount;
      my $D = $pij - $pi*$pj;
      $deltasum[$si-1] += ($D*$D)/($pi*(1-$pi)*($pj)*(1-$pj));
      $deltasum[-1] += ($D*$D)/($pi*(1-$pi)*($pj)*(1-$pj));
    }
    
  }
  if ((scalar(@out)) > 1000) {
  
    my $result = join "\n", @out;
    
    @out = ();
    
    unless ( open(SUMM, ">>$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit;
    }
    print SUMM "\n$result";
    close (SUMM);
    
  }
}

close (IN);

for (my $d = 0; $d < $retainedcount; $d++) {
  my $deltasum = sprintf "%.4f", $deltasum[$d];
  push @out, "$chrom\t$sites[$d]\t$rs[$d]\t$deltasum\t$buddies[$d]";
}

if ((scalar(@out)) > 0) {

  my $result = join "\n", @out;
  
  @out = ();
  
  unless ( open(SUMM, ">>$outfile") ) {
      print "Cannot open file \"$outfile\" to write to!!\n\n";
      exit;
  }
  print SUMM "\n$result";
  close (SUMM);

}

