#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_c $opt_f $opt_s $opt_v $opt_b $opt_w);

# Usage
my $usage = "
HetDifsPerSite.pl - get heterozygosity differences per site

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

Usage: perl HetDifsPerSite.pl options
  -c  chrom to examine
  -v  (path to) vcf (can be .gz)
  -b  (path to) outfile baseline name used in HetDifsPerRegion.pl (the only files starting with this name in the directory should be those you intend to use)
 optional:
  -f  min freq to include [default = 1%]
  -w  potential distance of a selection signal [default = 1000000]
  -s  subset of files to examine [site-site]
";

#############

getopts('c:o:f:s:v:b:w:');

die $usage unless ($opt_c);
die $usage unless ($opt_v);
die $usage unless ($opt_b);

my ($chrom, $minfreq, $genotypes, @subset, $baselinename, $relevantdist);

$chrom = $opt_c;

$genotypes = $opt_v;

$baselinename = $opt_b;

if (defined $opt_f) {
  $minfreq = $opt_f;
} else {
  $minfreq = 0.01;
}

if (defined $opt_w) {
  $relevantdist = $opt_w;
} else {
  $relevantdist = 1000000;
}

if (defined $opt_s) {
  @subset = split "-", $opt_s;
}

my @hetfiles = glob("$baselinename"."_$chrom"."_*");

unless (defined $hetfiles[0]) {
  print "No files to examine for $chrom!\n";
  exit;
}

my %pops;

my %HoHhets;

foreach my $hetfile (@hetfiles) {

  open(HET, $hetfile) || die "can't open $hetfile";
  
  print "Examining $hetfile.\n";
  
  my @region1 = split "$baselinename"."_$chrom"."_", $hetfile;
  
  my @region2 = split /\./, $region1[1];
  
  while (<HET>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[2] =~ /\d/) {
      my @ind1 = split "_", $data[0];
      my @ind2 = split "_", $data[1];
      $pops{$ind1[0]} = 1;
      $pops{$ind2[0]} = 1;
      my $combo = "$data[0]\t$data[1]";
      $HoHhets{$region2[0]}{$combo} = $data[2];
    }
  }
  
  close (HET);

}

my @pop = sort by_number (keys %pops);

my $mincount = $minfreq*2*(scalar(@pop));

print "Min count defined as $mincount.\n";

my @allsamples;

if ($genotypes =~ /.gz$/) {
  open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
}
else {
  open(IN, $genotypes) || die "can't open $genotypes";
}

my %HoHhetspergeno;

my %HoHcountspergeno;

my %freqs;

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
  if ((defined $subset[0])&&(($data[1] < $subset[0])||($data[1] > $subset[1]))) {
    next;
  }
  unless (defined $chrom) {
    $chrom = $data[0];
  }
  my @inds = @data[@allsamples[@pop]];
  my $flag = 0;
  my %seen;
  $seen{0} = 0;
  $seen{1} = 0;
  
  my $seenhemi = 0;
  if ($chrom =~ /\d/) {
    foreach my $i (@inds) {
      unless ($i =~ /\d\|\d/) {
        $flag = 1;
        last;
      }
      my @idata = split /\|/, $i;
      foreach my $id (@idata) {
        $seen{$id} +=1;
      }
    }
  } else {
    foreach my $i (@inds) {
      unless (($i =~ /^\d$/)||($i =~/\d\|\d/)) {
        $flag = 1;
        last;
      }
      if ($i =~ /^\d$/) {
        $seenhemi = 1;
        $seen{$i} +=1;
      } else {
        my @idata = split /\|/, $i;
        foreach my $id (@idata) {
          $seen{$id} +=1;
        }
      }
    }
  }
  if (($chrom !~ /\d/)&&($seenhemi == 0)) {
    $flag = 1;
  }
  if (($flag == 0)&&($seen{0} >= $mincount)&&($seen{1} >= $mincount)) {
    my $freq = sprintf "%.4f", $seen{1}/($seen{0}+$seen{1});
    if ($freq > 0.5) {
      $freq = 1 - $freq;
    }
    my $regiontouse;
    my $backupregion;
    foreach my $region (keys %HoHhets) {
      my @range = split "-", $region;
      if (($data[1] >= $range[0] + $relevantdist)&&($data[1] <= $range[1] - $relevantdist)) {
        $regiontouse = $region;
      } elsif (($data[1] >= $range[0])&&($data[1] <= $range[1])) {
        $backupregion = $region;
      } 
    }
    unless (defined $regiontouse) {
      if (defined $backupregion) {
        $regiontouse = $backupregion;
      }
    }
    if (defined $regiontouse) {
      my $site = "$data[1]\t$data[2]";
      $freqs{$site} = $freq;
      my @seeninds;
      my $ind1count = 0;
      foreach my $i1 (@inds) {
        my @i1data = split /\|/, $i1;
        my $ind2count = 0;
        if (defined $i1data[1]) {
          my $comboself = "$pop[$ind1count]"."_1\t$pop[$ind1count]"."_2";
          my @genoself = sort by_number ($i1data[0],$i1data[1]);
          my $genoself = join "\t", @genoself;
          if (defined $HoHhetspergeno{$site}{$genoself}) {
            $HoHhetspergeno{$site}{$genoself} += $HoHhets{$regiontouse}{$comboself};
            $HoHcountspergeno{$site}{$genoself} += 1;
          } else {
            $HoHhetspergeno{$site}{$genoself} = $HoHhets{$regiontouse}{$comboself};
            $HoHcountspergeno{$site}{$genoself} = 1;
          }
        }
        foreach my $i2 (@seeninds) {
          my @i2data = split /\|/, $i2;
          my $combo11 = "$pop[$ind2count]"."_1\t$pop[$ind1count]"."_1";
          my @geno11= sort by_number ($i2data[0],$i1data[0]);
          my $geno11 = join "\t", @geno11;
          if (defined $HoHhetspergeno{$site}{$geno11}) {
            $HoHhetspergeno{$site}{$geno11} += $HoHhets{$regiontouse}{$combo11};
            $HoHcountspergeno{$site}{$geno11} += 1;
          } else {
            $HoHhetspergeno{$site}{$geno11} = $HoHhets{$regiontouse}{$combo11};
            $HoHcountspergeno{$site}{$geno11} = 1;
          }
          if (defined $i1data[1]) {
            my $combo12 = "$pop[$ind2count]"."_1\t$pop[$ind1count]"."_2";
            my @geno12= sort by_number ($i2data[0],$i1data[1]);
            my $geno12 = join "\t", @geno12;
            if (defined $HoHhetspergeno{$site}{$geno12}) {
              $HoHhetspergeno{$site}{$geno12} += $HoHhets{$regiontouse}{$combo12};
              $HoHcountspergeno{$site}{$geno12} += 1;
            } else {
              $HoHhetspergeno{$site}{$geno12} = $HoHhets{$regiontouse}{$combo12};
              $HoHcountspergeno{$site}{$geno12} = 1;
            }
            if (defined $i2data[1]) {
              my $combo22 = "$pop[$ind2count]"."_2\t$pop[$ind1count]"."_2";
              my @geno22= sort by_number ($i2data[1],$i1data[1]);
              my $geno22 = join "\t", @geno22;
              if (defined $HoHhetspergeno{$site}{$geno22}) {
                $HoHhetspergeno{$site}{$geno22} += $HoHhets{$regiontouse}{$combo22};
                $HoHcountspergeno{$site}{$geno22} += 1;
              } else {
                $HoHhetspergeno{$site}{$geno22} = $HoHhets{$regiontouse}{$combo22};
                $HoHcountspergeno{$site}{$geno22} = 1;
              }
            }
          }
          if (defined $i2data[1]) {
            my $combo21 = "$pop[$ind2count]"."_2\t$pop[$ind1count]"."_1";
            my @geno21= sort by_number ($i2data[1],$i1data[0]);
            my $geno21 = join "\t", @geno21;
            if (defined $HoHhetspergeno{$site}{$geno21}) {
              $HoHhetspergeno{$site}{$geno21} += $HoHhets{$regiontouse}{$combo21};
              $HoHcountspergeno{$site}{$geno21} += 1;
            } else {
              $HoHhetspergeno{$site}{$geno21} = $HoHhets{$regiontouse}{$combo21};
              $HoHcountspergeno{$site}{$geno21} = 1;
            }
          }
          $ind2count +=1;
        }
        push @seeninds, $i1;
        $ind1count +=1;
      }
    }
  }
}

close (IN);

my @out;

foreach my $site (keys %HoHhetspergeno) {
  my $meanhet00;
  my $count00;
  if (defined $HoHhetspergeno{$site}{"0\t0"}) {
    $meanhet00 = sprintf "%.4f", $HoHhetspergeno{$site}{"0\t0"}/$HoHcountspergeno{$site}{"0\t0"};
    $count00 = $HoHcountspergeno{$site}{"0\t0"};
  }
  my $meanhet01;
  my $count01;
  if (defined $HoHhetspergeno{$site}{"0\t1"}) {
    $meanhet01 = sprintf "%.4f", $HoHhetspergeno{$site}{"0\t1"}/$HoHcountspergeno{$site}{"0\t1"};
    $count01 = $HoHcountspergeno{$site}{"0\t1"};
  }
  my $meanhet11;
  my $count11;
  if (defined $HoHhetspergeno{$site}{"1\t1"}) {
    $meanhet11 = sprintf "%.4f", $HoHhetspergeno{$site}{"1\t1"}/$HoHcountspergeno{$site}{"1\t1"};
    $count11 = $HoHcountspergeno{$site}{"1\t1"};
  }
  if ((defined $meanhet00)&&(defined $meanhet11)) {
    my $dif = abs($meanhet00-$meanhet11);
    push @out, "$site\t$freqs{$site}\t$meanhet00\t$count00\t$meanhet01\t$count01\t$meanhet11\t$count11\t$dif";
  }
}

my $result = join "\n", @out;

@out = ();

my $outfile = "$baselinename"."_$chrom"."_Difs_Per_Site.txt";

unless ( open(SUMM, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print SUMM "Site\tID\tFreq\tHet00\tCount00\tHet01\tCount01\tHet11\tCount11\tDif\n$result";
close (SUMM);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
