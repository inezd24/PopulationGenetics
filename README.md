## IBD analyses
Last updated: Tuesday 31 October 2023
In this script, you will PLINK (1.9) data to phase with BEAGLE, and subsequently use for hap-ibd and IBDne analyses.
The script can be run in bash and plots can be recreated in R with the ggplot2 package. 

### Phasing with BEAGLE
Script: beagle_raute_only.sh or beagle_data.sh. See also beagle_example.sh from developers.

1. Convert files to VCF format.
```
plink --bfile ~/plink/our_samples/our_samples_maf --allow-no-sex --recode vcf-iid --alleleACGT --out ~/plink/phased_data/our_samples
```

2. Compress the VCF files and index it.
```
bgzip ~/plink/phased_data/our_samples.vcf
tabix -f -p vcf ~/plink/phased_data/our_samples.vcf.gz
```

3. Create beagle directory, move into directory, and download beagle. Also create new directory for specific files.
```
mkdir ./beagle
cd ~/beagle
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
mkdir raute_only
```

4. Download GRCh37 maps in PLINK format from Browning lab, unzip, and rename X chromosome to 23.
```
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
unzip plink.GRCh37.map.zip
mv plink.chrX.GRCh37.map plink.chr23.GRCh37.map
````

5. Split vcf file by chromosomes (see Epi's repository for source code: https://github.com/epifaniarango/popgen_with_epi/tree/IBD).
```
for chr in `bcftools view -h ~/plink/phased_data/our_samples.vcf.gz | perl -ne '
 if (/^##contig=<ID=([^,]+)/) { if (length($1)<=2) {
   print "$1\n"
 } }'`; do
  bcftools view -Oz -r $chr ~/plink/phased_data/our_samples.vcf.gz  > ./raute_only/our_samples$chr.vcf.gz &
done
```

6. Run BEAGLE for chromosomes separately.
```
for chromosome in {1..22}; do
      java -jar beagle.22Jul22.46e.jar gt=./raute_only/our_samples${chromosome}.vcf.gz  map=plink.chr${chromosome}.GRCh37.map         window=20 out=./raute_only/our_samples_chr${chromosome} nthreads=12
done
```

7. Index each chromosome vcf for later analyses.
```
for chromosome in {1..22}; do
	bgzip ./raute_only/our_samples_chr${chromosome}.vcf
	tabix -f -p vcf ./raute_only/our_samples_chr${chromosome}.vcf.gz
done
```

8. Concatenate the chromosome-specific VCF files, compress, and index (for downstream analyses).
```
bcftools concat -o concat_raute.vcf our_samples_chr{1..22}.vcf
bgzip concat_raute.vcf
tabix -p vcf concat_raute.vcf.gz
```

9. Lastly, also merge the genetic maps together (still in PLINK format) for downstream analyses.
```
cat plink.chr{1..22}.GRCh37.map > GRCh37.map
```

### Identifying IBD (identity-by-descent) segments using hap-IBD
Script: in IBD_analyses directory, find script: hap_IBD.sh. See also run_hap_IBD_test.sh and test folder for examples. 

1. Create directory for all IBD analyses and download hap-ibd program.
```
cd; mkdir IBD_analyses; cd ~/IBD_analyses
wget https://faculty.washington.edu/browning/hap-ibd.jar
```

2. Create directory for hap-ibd output.
```
mkdir hap_ibd
mkdir hap_ibd/single_run #as I previously ran BEAGLE multiple times
```

3a. Run hap-ibd separated for each chromosome
```
for chromosome in {1..22};
do java -jar ./hap-ibd.jar  gt=~/plink/phased_data/phased_no_LD/single_run/nonLD_phased${chromosome}.vcf.gz map=~/beagle/plink.chr${chromosome}.GRCh37.map out=./hap_ibd/single_run/hap_ibd${chromosome}.out
done
```

3b. Run hap-ibd with the combined VCFs (the concatenated file)
```
java -jar ./hap-ibd.jar gt=~/beagle/raute_only/concat_raute.vcf.gz map=~/beagle/GRCh37.map out=./hap_ibd/concat_raute

```

### Estimate effective population sizes with IBDne.
Script: in IBD_analyses directory, find scripts: extract_raute_IBD.sh and IBD_ne.sh.

1. Download IBDne program in the previously made directory.
```
cd ~/IBD_analyses
wget https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
```

2. Unzip the hap-ibd output.
```
gunzip ~/IBD_analyses/hap_ibd/single_run/*.ibd.gz 
gunzip ~/IBD_analyses/hap_ibd/concat_raute.ibd.gz
```

3. Since the hap-ibd software ran for all our samples, 1342 in total at this time (31 oct 2023), we need to extract only the Raute samples.

3a. Select all individuals from the full census list (of our study population only) and select only the names of the individuals to be put in a new file (Raute_individuals.txt).
```
grep "Raute" ~/plink/modified_samples/all_individuals.txt > 'IBD_Raute.txt'
awk '{print $2}' IBD_Raute.txt > Raute_individuals.txt
```

3b. Take the selection previously made and apply it to the .ibd output files from hap-ibd, so that only the segments between Raute individuals remain. Remove remaining unnecessary files after inspection.
```
for chr in {1..22}; do
	awk 'NR==FNR{a[$0];next}($1 in a)' Raute_individuals.txt hap_ibd/single_run/hap_ibd${chr}.out.ibd > hap_ibd/single_run_raute/raute_ibd${chr}.out
	awk 'NR==FNR{a[$0];next}($3 in a)' Raute_individuals.txt hap_ibd/single_run_raute/raute_ibd${chr}.out > hap_ibd/single_run_raute/raute_ibd${chr}
	rm hap_ibd/single_run_raute/raute_ibd${chr}.out
done
```

3c. Do the same for the concatenated file.
```
awk 'NR==FNR{a[$0];next}($1 in a)' Raute_individuals.txt ./hap_ibd/concat_raute.ibd > ./hap_ibd/concat_raute.ibd2
awk 'NR==FNR{a[$0];next}($1 in a)' Raute_individuals.txt ./hap_ibd/concat_raute.ibd2 > ./hap_ibd/concat_raute_113.ibd
rm ./hap_ibd/concat_raute.ibd2
```

4a. Run IBDne on chromosome specific files. Note: you are required to put 'nboots=0', because bootstrapping is not possible when you run IBDne per chromosome. Hence, this will give you some estimates of how the results might look, but it will not give you confidence intervals. Here, I also didn't specify a minimum number of centimorgans that the segments should have in order to be considered. 
```
for chr in {1..22}; do
cat ~/IBD_analyses/hap_ibd/single_run_raute/raute_ibd${chr} | java -jar ~/IBD_analyses/ibdne.23Apr20.ae9.jar map=~/beagle/plink.chr${chr}.GRCh37.map out=./ne_raute/ne_${chr}.out nthreads=12 nboots=0
done
```

4b. Run IBDne on the concatenated file. This time, I use nboots=100 (instead of the default 80) and mincm=7 (so minimum of 7 centimorgans). You can also run it once without specifying mincm (the default is 2) but with bootstrapping to see how changing the cM threshold changes the results. 
```
cat ~/IBD_analyses/hap_ibd/concat_raute_113.ibd | java -jar ~/IBD_analyses/ibdne.23Apr20.ae9.jar map=~/beagle/GRCh37.map out=./ne_raute/concat_all nthreads=12 nboots=100 mincm=7
```

### Visualisation in R
Now you can analyse the results in R using a simple ggplot plot and some small modifications:
```
# Load (and if necessary download) packages
library(ggplot2)

# Read file
concat.out <- read.delim("~/PHD/Genetic_Data/IBD_analyses/raute_all_concat/concat_all.ne")

# log10 transform the ne estimates, lower CI and upper CI. 
concat.out$log10_ne <- log10(concat.out$NE)
concat.out$log10_lwr <- log10(concat.out$LWR.95.CI)
concat.out$log10_upr <- log10(concat.out$UPR.95.CI)

# Select only the last 50 generations, since this is where the program is most reliable. 
concat.out <- concat.out[!(concat.out$GEN >50),]

# Plot ne for these 50 generations
ggplot(concat.out, aes(x=GEN, y=log10_ne)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log10_lwr, ymax = log10_upr), 
              alpha=0.1, 
              linetype="dashed",
              color="grey") +
  scale_color_grey() + 
  theme_classic()
```


