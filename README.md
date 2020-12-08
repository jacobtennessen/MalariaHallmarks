# MalariaHallmarks
Scripts for detecting signatures of natural selection characteristic of human genes that impact resistance to malaria.

These scrips start with phased genotype data in vcf format as in the 1000 Genomes data, with a separate file for each chromosome. Other variations of the vcf format may not work.

To calculate Dng (aka 🍡) for a test of Divergent Haplotypes:

1.	Use GetDango.pl on vcf data.

e.g. perl GetDango.pl -v ALL.chr1_GRCh38.genotypes.vcf.gz -a 0,3,4,7,8,9 -o Dango_chr1.txt -w 5000

#Do this for every vcf file.

To calculate TR (aka 🔱) for a test of Repeated Shifts:

1.	Use FstPerSite.pl on vcf data to calculate pairwise FST values between populations.

e.g. perl FstPerSite.pl -m /path/to/directory/ -p 0,3,4,7,8,9-1,2,5,6,10 -v ALL.chrCHROM_GRCh38.genotypes.vcf.gz -f 0.01 -o Fst_values_over_0.01_pops_1_2.txt

#Do this three times, for all three pairs of three populations. All vcf files can be analyzed in a single job if they have the same name except for the chromosome number, which is indicated with the CHROM placeholder in the name.

2.	Use RegionsFstForTrident.pl on the output of FstForTrident.pl.

e.g. perl RegionsFstForTrident.pl -a Fst_values_over_0.01_pops_1_2.txt -r 50000 -o Fst_Regions_50k_pops_1_2.txt

#Do this three times, for all three pairs of three populations.

3.	Use FindViableFstRegionsForTrident.pl on the vcf data to find how many genomic windows could potentially be outliers

e.g. perl FindViableFstRegionsForTrident.pl -m /path/to/directory/ -a 0,3,4,7,8,9 -b 1,2,5,6,10 -c 11,12,13,14,15 -v ALL.chrCHROM_GRCh38.genotypes.vcf.gz -r 50000 -o Viable_Fst_Regions_50k.txt

#Again, use the generic vcf name with CHROM placeholder.

4.	Use GetTrident.pl to calculate the final statistic.

e.g. perl GetTrident.pl -a Fst_Regions_50k_pops_1_2.txt -b Fst_Regions_50k_pops_2_3.txt -c Fst_Regions_50k_pops_1_3.txt -r 50000 -v Viable_Fst_Regions_50k.txt -o Trident_50k_pops_1_2_3.txt

To calculate ΠAHz (aka ⏸️) for a test of Arrested Sweep:

1.	Use HetDifsPerRegion.pl on vcf data to calculate heterozygosity differences per pair of individuals

e.g perl HetDifsPerRegion.pl -m /path/to/directory/ -a 0,3,4,7,8,9 -v ALL.chrCHROM_GRCh38.genotypes.vcf.gz -d Chrom_Data.txt -c 2 -o Het_Difs_Pop1

#For efficiency, this can optionally by run on each vcf separately, or sets of vcfs, in parallel (e.g. if each chromosome is a separate vcf). As with TR/🔱, use the generic vcf name with CHROM placeholder.

2.	Use HetDifsPerSite.pl on vcf data to calculate average heterozygosity differences for each site.

e.g. perl HetDifsPerSite.pl -c 2 -v /path/to/directory/ALL.chr2_GRCh38.genotypes.vcf.gz -b Het_Difs_Pop1

#This must be run vcf-by-vcf. So unlike some of the other scripts, one enters the actual vcf name together with its path. Do this for each vcf.

3.	Use FstPerSite.pl on vcf data to calculate FST values among ingroup populations.

e.g. perl FstPerSite.pl -m /path/to/directory/ -v ALL.chrCHROM_GRCh38.genotypes.vcf.gz -p 0,3,4,7,8,9-1,2,5,6,10-11,12,13,14,15-16,19,21,22,23,25-17,18,20,24,26 -o Fst_Within_Ingroup.txt

#Use the generic vcf name with CHROM placeholder

4.	Use FstMatrix.pl on the output of FstPerSite.pl

e.g. perl FstMatrix.pl -f Fst_Within_Ingroup.txt -x
-o Fst_Within_Ingroup_Autosomes_Matrix.txt 

5.	Use FreqPerSite.pl on vcf data to get outgroup frequencies.

e.g. perl FreqPerSite.pl -m /path/to/directory/ -p 27,28,29,30,31,32 -v ALL.chrCHROM_GRCh38.genotypes.vcf.gz -o Freqs_outgroup.txt

#Use the generic vcf name with CHROM placeholder.

6.	Use GetPause.pl with the output of HetDifsPerSite.pl and FstMatrix.pl

e.g. perl GetPause.pl -b Het_Difs_Pop1 -f Fst_Within_Ingroup.txt -i Fst_Within_Ingroup_Autosomes_Matrix.txt -d Chrom_Data.txt -e Freqs_outgroup.txt -m 10000 -x -o Pause_output.txt

