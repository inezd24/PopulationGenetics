# Finestructure 
Pipeline to use the program finestructure with merged genome-wide SNP array data from differented populations. This requires phasing data, which you can do with several programs. I have used BEAGLE (per recommendation of others + ease of use) and you can find the pipeline to do this under the branch 'IBD_analyses'. This software does not like missing data, hence filter for missing genotypes (and SNPs) as well as minor allele frequencies prior to using your data. 

### Preparing files for chromopainter input (test run)
Necessary scripts for this part are try_out_chrom.sh and chr_EM10.sh as well as the perl conversion scripts vcf2cp.pl and convertrecfile.pl. On my server, a lot of the downloading I have not put in a script but done separately. I've tried to add most of it here as otherwise I won't be able to replicate it, but some may be missing (e.g. the prereqs for the GUI took me quite a while to install). 

1. Download, rename, and install finestructure
```
cd; wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.1.1.zip
unzip fs_4.1.1.zip; rm fs_4.1.1.zip
mv fs_4.1.1 finestructure
cd finestructure
./fs_install.sh
```

2. OPTIONAL: You can download the GUI, although it is depracated. This means you should run analyses from the pipeline, but may use the GUI (unsupported). Just in case, this is how to install the GUI. It has several dependencies: GNU Scientific Library (GSL, development version), automake (package automake in Ubuntu and OpenSUSE). wxwidgets 3.0.0+ (packages libwxgtk3.0-dev and wx-common in Ubuntu, wxGTK-devel in OpenSUSE; WARNING: This has many dependencies.), and GCC c++ compiler. If you need to install any of these things, check the respective packages and programs. Also upgrade packages you already have on ubuntu to make sure everything has the latest version and Ubuntu has the latest version as well. I don't use the GUI after this, hence this step is optional and just for exploration. 
```
wget https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure-0.1.0GUI.tar.gz
tar -xzvf finestructure-0.1.0GUI.tar.gz # (check correct version)
cd finestructure-0.1.0/gui
./configure
make
sudo make install # optional; installs finegui exceutable
```

3. Download and install chromopainterv2
```
wget https://github.com/sahwa/ChromoPainterV2/blob/main/ChromoPainterv2
wget https://github.com/sahwa/ChromoPainterV2/blob/main/ChromoPainterv2.c
gcc -Wall -o ChromoPainterv2 ChromoPainterv2.c -lm -lz
rm ChromoPainterv2.c
```

4. Download and unzip references files/genetic maps. 
```
cd; git clone https://github.com/adimitromanolakis/geneticMap-GRCh37
gunzip ~/geneticMap-GRCh37/genetic_map_GRCh37_*
```

5. Make test directory and enter directory. 
```
cd finestructure; mkdir chr1; cd ./chr1
```

6. Copy BEAGLE phased files into this working directory
```
for chr in {1,5,8,15}; do
cp ~/plink/phased_data/phased_no_LD/single_run/nonLD_phased${chr}.vcf.gz ./
gunzip nonLD_phased${chr}.vcf.gz
done
```

7. Convert these files to chromopainter format via the safe VCF raute provided by finestructure (vcf2cp.pl). See https://www.biostars.org/p/9519621/ or https://github.com/danjlawson/finestructure4 for extra help and the necessary code. 
```
for chr in {1,5,8,15}; do
perl ~/finestructure/vcf2cp.pl nonLD_phased${chr}.vcf nonLD_phased${chr}
done
```

8. Make recombination file with convertrecfile.pl. Do this for a subset of chromosomes only.
```
for chr in {1,5,8,15}; do
~/finestructure/convertrecfile.pl -M hap ~/finestructure/chr1/nonLD_phased${chr}.phase ~/geneticMap-GRCh37/genetic_map_GRCh37_chr${chr}.txt chr${chr}.recombfile 
done
```

9. Run a test with chromopainterv2 for the selected chromosomes. This will be necessary for when we run chromopainter for all chromosomes. 
```
for chr in {1,5,8,15}; do
~/finestructure/ChromoPainterv2 -g nonLD_phased${chr}.phase -r chr${chr}.recombfile -t nonLD_phased${chr}.ids -o chr${chr}_EM10 -a 0 0 -i 10 -in -iM
done
```

10. Using the '.EMprobs.out' file, calculate the average switch parameter and (global) mutation rate in R.
```
# Script executable in R
# Developed by Inez Derkx
# Created on Friday 3 November 2023
# Last updated on Monday 6 November 2023

################################################# START OF SCRIPT #################################################

# Import libraries
library(dplyr)

# Import files
chr5_EM10.EMprobs <- read_table("~/PHD/Genetic_Data/finestructure/chr5_EM10.EMprobs.out_notitle", col_names = F)
chr8_EM10.EMprobs <- read_table("~/PHD/Genetic_Data/finestructure/chr5_EM10.EMprobs.out_notitle", col_names = F)
chr15_EM10.EMprobs <- read_table("~/PHD/Genetic_Data/finestructure/chr5_EM10.EMprobs.out_notitle", col_names = F)

# Function to change files
wrangle_EM_probs <- function(x, nr){
  
  # first delete lines we dont want
  x <- x[!is.na(x$X2),]
  
  # Then add a sixth column to identify which chromosome
  x$chr <- rep(nr, each = nrow(x))
  
  # Then add id's of the individuals (not actual ids, just replacements)
  x <- x %>%
  mutate(ind = cumsum(X1 == 0))
  
  # Now rename the first 5 columns
  names(x)[1:5] <- c("EM", "loglike_snp_hap1", "loglike_snp_hap2", "switch_rate", "global_mut_rate")
  
  # Return
  return(x)
  
  }

# Apply function
chr5_EM10.EMprobs <- wrangle_EM_probs(chr5_EM10.EMprobs, 5)
chr8_EM10.EMprobs <- wrangle_EM_probs(chr8_EM10.EMprobs, 8)
chr15_EM10.EMprobs <- wrangle_EM_probs(chr15_EM10.EMprobs, 15)

# Concatenate files
EM10.EMprobs <- rbind(chr5_EM10.EMprobs, chr8_EM10.EMprobs, chr15_EM10.EMprobs)

# For mutation probability
EM10_mut <-EM10.EMprobs %>% group_by(ind, chr) %>% summarise(mut_rate = mean(global_mut_rate))
EM10_mut <-EM10_mut %>% group_by(ind) %>% summarise(mut_rate = mean(mut_rate))
mean(EM10_mut$mut_rate) #0.0007994967

# For switch rate
EM10_switch <-EM10.EMprobs %>% group_by(ind, chr) %>% summarise(switch_rate = mean(switch_rate))
EM10_switch <-EM10_switch %>% group_by(ind) %>% summarise(switch_rate = mean(switch_rate))
mean(EM10_switch$switch_rate) #123.0597

################################################## END OF SCRIPT ##################################################
```

11. Now that we have run a test, we can run chromopainter for all chromosomes. For this, we first have to copy the previous steps 5-8 for all chromosomes. Here, they are put together in one script:
```
# Make new directory for new input files and output files
mkdir ~/finestructure/output_chromo
mkdir ~/finestructure/all_chrs

# Copy files from previous run into new directory (we don't have to make new recombination/chromopainter input files for these of course). 
cp ~/finestructure/chr1/nonLD_phased* ~/finestructure/all_chrs/
cp ~/finestructure/chr1/*.recombfile ~/finestructure/all_chrs/

# Copy remaining necessary vcfs to the new directory
for chr in {2,3,4,6,7,9,10,11,12,13,14,16,17,18,19,20,21,22}; do
	cp ~/plink/phased_data/phased_no_LD/single_run/nonLD_phased${chr}.vcf.gz ~/finestructure/all_chrs/
	gunzip ~/finestructure/all_chrs/nonLD_phased${chr}.vcf.gz
done

# Convert it to chromopainter format via the safe VCF route:
for chr in {2,3,4,6,7,9,10,11,12,13,14,16,17,18,19,20,21,22}; do
	perl ~/finestructure/vcf2cp.pl ~/finestructure/all_chrs/nonLD_phased${chr}.vcf ~/finestructure/all_chrs/nonLD_phased${chr}
done

# Make recombination map with convertrecfile.pl
for chr in {2,3,4,6,7,9,10,11,12,13,14,16,17,18,19,20,21,22}; do
	~/finestructure/convertrecfile.pl -M hap ~/finestructure/all_chrs/nonLD_phased${chr}.phase ~/geneticMap-GRCh37/genetic_map_GRCh37_chr${chr}.txt ~/finestructure/all_chrs/chr${chr}.recombfile
done
```

12. Run chromopainterv2 with found parameters and for all chromosomes. The switchrate is given with -n and the global mutation rate is given with -M. We use '-a 0 0' to indicate we want to paint all individuals. Usually we'd have to specify an -f file, but not if you use the -a option. The values used here for -M and -n are obtained through the R script above. 
```
for chr in {1..22}; do
	~/finestructure/ChromoPainterv2 -g ~/finestructure/all_chrs/nonLD_phased${chr}.phase -r ~/finestructure/all_chrs/chr${chr}.recombfile -t ~/finestructure/all_chrs/nonLD_phased${chr}.ids -o ~/finestructure/output_chromo/chr${chr}_final -a 0 0 -n 123.0597 -M 0.0007994967
done
```






