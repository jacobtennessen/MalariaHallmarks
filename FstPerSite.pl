#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_m $opt_v $opt_o $opt_p $opt_f $opt_q $opt_t $opt_r $opt_l $opt_s);

# Usage
my $usage = "
FstPerSite.pl - Calculates Fst from original vcfs

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

Usage: perl FstPerSite.pl options
  -m  metafolder
  -v  generic vcf or vcf.gz name with CHROM, e.g. ALL.chrCHROM_GRCh38.genotypes.vcf.gz
                (assumes each chromosome is a separate vcf in metafolder,
                all with the same name expect for the chromosome,
                e.g. ALL.chr1_GRCh38.genotypes.vcf.gz, ALL.chr2_GRCh38.genotypes.vcf.gz, etc.)
  -p  '-'-delimited list of comma-delimited lists of positions (samples start at 0) (e.g. 0,3,4,7,8,9-1,2,5,6,10 for two populations of 6 and 5 individuals)
  -o  (path to) outfile for list of high-Fst SNPs
 optional:
  -f  minimum Fst for export to candidate list [default = none]
  -t  maximum Fst for export to canddiate list [default = none]
  -q  minimum MAF for export to candidate list [default = 0]
  -r  maximum MAF for export to candidate list [default = 0.5]
  -l  comma-delimited list of chromosomes to examine (default = X,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
  -s  number of lines to output at once [default = entire chromosome]
";

#############

getopts('m:v:o:p:f:t:q:r:l:s:');

die $usage unless ($opt_m);
die $usage unless ($opt_o);
die $usage unless ($opt_p);
die $usage unless ($opt_v);

my ($metafolder, $outfile, @allpops, $fstthresh, $fstmax, $minfreq, $maxfreq, $vcfname, $outputsize, @chroms);

$metafolder = $opt_m;

$outfile = $opt_o;

@allpops = split "-", $opt_p;

$vcfname = $opt_v;

if (defined $opt_l) {
  @chroms = split ",", $opt_l;
} else {
  @chroms = ("X",1..22);
}

my %HoApops;

my $popcount = 0;

foreach my $plist (@allpops) {
  $popcount +=1;
  @{$HoApops{$popcount}} = split ",", $plist;
}

if (defined $opt_f) {
  $fstthresh = $opt_f;
}

if (defined $opt_t) {
  $fstmax = $opt_t;
}

if (defined $opt_q) {
  $minfreq = $opt_q;
} else {
  $minfreq = 0;
}

if (defined $opt_r) {
  $maxfreq = $opt_r;
} else {
  $maxfreq = 0.5;
}

if (defined $opt_s) {
  $outputsize = $opt_s;
}

my @allsamples;

my $seensnps = 0;

my @titlefreqs;

for (my $p = 1; $p <= $popcount; $p++) {
  push @titlefreqs, "Freq$p";
}

my $titlefreqs = join "\t", @titlefreqs;

my @titlesizes;

for (my $p = 1; $p <= $popcount; $p++) {
  push @titlesizes, "Size$p";
}

my $titlesizes = join "\t", @titlesizes;

unless ( open(OFC, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OFC "Chrom\tSite\tSNP\tRef\tAlt\t$titlefreqs\t$titlesizes\tFreq\tFst";
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
    my %HoAinds;
    for (my $p = 1; $p <= $popcount; $p++) {
      @{$HoAinds{$p}} = @data[@allsamples[@{$HoApops{$p}}]];
    }
    my $genoflag = 0;
    my %seenref;
    my %seenalt;
    my $seenref = 0;
    my $seenalt = 0;
    for (my $pop = 1; $pop <= $popcount; $pop ++) {
      $seenref{$pop} = 0;
      $seenalt{$pop} = 0;
      foreach my $i1 (@{$HoAinds{$pop}}) {
        if ($i1 =~ /0\|0/) {
          $seenref{$pop} += 2;
          $seenref +=2;
        } elsif (($i1 =~ /0\|1/)||($i1 =~ /1\|0/)||($i1 =~ /2\|0/)||($i1 =~ /0\|2/)) {
          $seenref{$pop} += 1;
          $seenalt{$pop} += 1;
          $seenref += 1;
          $seenalt += 1;
        } elsif (($i1 =~ /1\|1/)||($i1 =~ /1\|2/)||($i1 =~ /2\|1/)||($i1 =~ /2\|2/)) {
          $seenalt{$pop} += 2;
          $seenalt += 2;
        } elsif ($i1 =~ /^0$/) {
          $seenref{$pop} +=1;
          $seenref +=1;
        } elsif (($i1 =~ /^1$/)||($i1 =~ /^2$/)) {
          $seenalt{$pop} +=1;
          $seenalt +=1;
        } else {
          $genoflag += 1;
        }
      }
    }
    if (($genoflag == 0)&&($seenref > 0)&&($seenalt > 0)) {
      my @sizes;
      my @freqs;
      my $totalsize = 0;
      my $totalcount = 0;
      my @roundfreqs;
      for (my $pop = 1; $pop <= $popcount; $pop ++) {
        my $size = ($seenalt{$pop}+$seenref{$pop});
        push @sizes, $size;
        $totalsize += $size;
        my $freq = $seenalt{$pop}/$size;
        push @freqs, $freq;
        $totalcount += $seenalt{$pop};
        push @roundfreqs, sprintf "%.4f", $freq;
      }
      my $meanfreq = sprintf "%.4f", $totalcount/$totalsize;
      if ($meanfreq > 0.5) {
        $meanfreq = 1 - $meanfreq;
      }
      if (($meanfreq >= $minfreq)&&($meanfreq <= $maxfreq)) {
        my $sizes = join ",", @sizes;
        my $freqs = join ",", @freqs;
        my $fst = sprintf "%.4f", (Fstmulti ($sizes, $freqs));
        unless ((defined $fstthresh)&&($fst < $fstthresh)) {
          unless ((defined $fstmax)&&($fst > $fstmax)) {
            my $freqlist = join "\t", @roundfreqs;
            my $sizelist = join "\t", @sizes;
            push @outliers, "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$freqlist\t$sizelist\t$meanfreq\t$fst";
          }
        }
      }
      $seensnps +=1;
    }
    if ((defined $outputsize)&&((scalar (@outliers)) >= $outputsize)) {
      my $outliers = join "\n", @outliers;
      @outliers = ();
      unless ( open(OFT, ">>$outfile") ) {
          print "Cannot open file \"$outfile\" to write to!!\n\n";
          exit;
      }
      print OFT "\n$outliers";
      close (OFT);
    }
  }
  
  close (IN);
  
  if ((scalar (@outliers)) >= 1) {
  
    my $outliers = join "\n", @outliers;
    
    @outliers = ();
    
    unless ( open(OFC, ">>$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit;
    }
    print OFC "\n$outliers";
    close (OFC);
  
  }
  
  print "Chrom $c complete; $seensnps SNPs seen.\n";

}

###################

sub Fstmulti {
    my ($sizes, $freqs) = @_;
    
    my @freqs = split ",", $freqs;
    my @samplesizes = split ",", $sizes;
    
    my %seenfreqs;
    
    foreach my $f (@freqs) {
      $seenfreqs{$f} = 1;
    }
    
    my $freqcount = scalar (keys %seenfreqs);
  
    my $Fst;
    
    if (($freqcount <= 1)&&((defined $seenfreqs{0})||(defined $seenfreqs{1}))) {
      $Fst = 0; 
    } else {
    
      my $num_sub_pops = scalar(@samplesizes);
  
      my ($TS_sub1,$TS_sub2);
      my @alleles = (0,1);
  
      my $avg_samp_size         = 0;
      my $avg_allele_freq       = 0;
      my $total_samples_squared = 0;
  
      foreach (my $c = 0; $c < $num_sub_pops; $c ++) {
        my $s = $samplesizes[$c];
        $avg_samp_size += $s;
        $total_samples_squared += $s**2;
        my $all_freq = $freqs[$c];
        $avg_allele_freq += $s * $all_freq;
      }
      
      my $total_samples =  $avg_samp_size;
      $avg_samp_size /= $num_sub_pops;
      $avg_allele_freq /= $total_samples;
  
      my $adj_samp_size = ( 1/ ($num_sub_pops - 1)) * ( $total_samples - ( $total_samples_squared/$total_samples));
  
      my $variance              = 0;
      my $sum_variance          = 0;
      my $i = 0;
      for (my $d = 0; $d < $num_sub_pops; $d ++) {
        my $s = $samplesizes[$d];
        $sum_variance += $s * (($freqs[$d] - $avg_allele_freq)**2);
      }
      $variance = ( 1 / (( $num_sub_pops-1)*$avg_samp_size))*$sum_variance;
  
      $TS_sub1 =
        $variance -
        (
          (1/($avg_samp_size-1))*
          (
            ($avg_allele_freq*(1-$avg_allele_freq)) -
            ( (($num_sub_pops-1)/$num_sub_pops)*$variance)
          )
        );
         
      $TS_sub2 =
        (
          (($adj_samp_size-1)/($avg_samp_size-1))*$avg_allele_freq*(1-$avg_allele_freq)
        ) +
        (
          1 +
          ( (($num_sub_pops-1)*($avg_samp_size-$adj_samp_size))/ ($avg_samp_size - 1))
        ) * 
        ($variance/$num_sub_pops);
  
          $Fst = $TS_sub1 / $TS_sub2;
    }
   
    return $Fst;
}

