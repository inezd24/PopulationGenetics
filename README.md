# PLINK to examine, filter and prune (LD) SNP array data
Here, I am using data genotyped on the Human Origins array (GWSA with 600.000 SNPs developed for population geneticicists). I am using PLINK version 1.9 (https://www.cog-genomics.org/plink/) unless specified differently (e.g., 'plink2'). We also use IBIS and R to extract unrelated individuals based on a kinship coefficient obtained through IBD segments.

### Initial filtering and pruning of samples
Use of script 'filter_our_samples.sh' from PLINK directory

1. Download and install PLINK. Make a new directory for storing the files used for/created with PLINK and move into directory. 
```
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip #Download 64-bit version for Linux
sudo unzip plink_linux_x86_64_20231018.zip -d plink_install # Unzip
cd plink_install #install
sudo cp plink /usr/local/bin # make accessible from command line
sudo chmod 755 /usr/local/bin/plink # make accessible from command line
sudo nano ~/.bashrc # add PLINK to path
# This will open the nano editor in the terminal. Scroll down to the bottom of the file and add the following line:
# export PATH=/usr/local/bin:$PATH
# Click cmd+s and cmd+x to save and exit the file.
# To test the installation, restart the terminal and run: plink
mkdir ~/plink; cd plink
```

2. Use PLINK 1.9 to filter VCF file for only the Raute samples and convert to plink format
```
## FILTER FOR RAUTE
VCF=VCF/GenotypingRaute2.vcf
Raute_INCL=~/plink/our_samples/update_data/samples_raute.txt
SEX=~/plink/our_samples/update_data/update_sex_ids_raute.txt # Here I add the sex ids based on what we know from the census data

# use PLINK 1.9 to filter the file and only keep Raute
plink --vcf $VCF --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-missing-var-ids @:# \
        --keep $Raute_INCL \
	--update-sex $SEX \
        --geno 0.1 \
        --mind 0.1 \
        --hwe 0.000001 midp \
        --make-bed \
        --out ./our_samples/ibis/raute_samples

## FILTER FOR DAILEKH SAMPLES
VCF=VCF/GenotypingRaute2.vcf
Controls_INCL=~/plink/our_samples/update_data/samples_controls.txt
SEX=~/plink/our_samples/update_data/update_sex_ids_controls.txt # Here I add the sex ids based on what we know from the census data

# use PLINK 1.9 to filter the file and only keep Raute
plink --vcf $VCF --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-missing-var-ids @:# \
        --keep $Controls_INCL \
	--update-sex $SEX \
        --geno 0.1 \
        --mind 0.1 \
        --hwe 0.000001 midp \
        --make-bed \
        --out ~/plink/our_samples/ibis/controls_samples

```

3. Use plink 2 to set variant ids as 'chr_position' (plink 1.9 doesn't do this) and remove duplicates
```
## FILTER FOR RAUTE
PHENO=~/plink/our_samples/update_data/update_raute_2023_pheno.txt # Here we update phenotypes so all have the same for downstream analyses
IDS=~/plink/our_samples/update_data/update_raute_2023_IDS.txt # here I update the names to 'Population<nr>', like for instance 'Agta1', 'Agta2', etc. and update population names. 
plink2 --bfile ./our_samples/ibis/raute_samples --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
        --pheno $PHENO \
        --update-ids $IDS \
	--make-bed \
	--out ./our_samples/ibis/raute_ibis

rm ./our_samples/ibis/raute_samples.*

## FILTER FOR DAILEKH SAMPLES
PHENO=~/plink/our_samples/update_data/update_controls_2023_pheno.txt
IDS=~/plink/our_samples/update_data/update_controls_2023_IDS.txt
plink2 --bfile ./our_samples/ibis/controls_samples --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
        --pheno $PHENO \
        --update-ids $IDS \
	--make-bed \
	--out ./our_samples/ibis/controls_ibis

rm ~/plink/our_samples/ibis/controls_samples.*
```


### Using IBIS to extract kinship coefficient and determine unrelated individuals in R
See https://github.com/williamslab/ibis/blob/master/README.md for more instructions on how to use IBIS as well as an overview of all the possible parameters. IBIS is aimed at unphased datasets. We use IBIS on the datafiles that are not yet filtered for MAF or pruned for LD. 

1. Compile IBIS
```
cd
git clone --recurse-submodules https://github.com/williamslab/ibis.git
cd ibis; make
git pull
git submodule update --remote # to pull updates
```

2. Run IBIS
```
cd ~/plink/our_samples/ibis/

# Locate files to use
find ~/plink/our_samples/ibis/ -type f -name '*.bed' | awk -F/ '{print $NF}' | cut -d. -f1 | sort -u > populations

# Perform IBIS for each population
while IFS= read -r file; do
	~/ibis/add-map-plink.pl ${file}.bim ~/ibis/genetic_map_GRCh37_chr{1..22}.txt > ~/plink/our_samples/ibis/${file}_map.bim
	~/ibis/ibis ${file}.bed ${file}_map.bim ${file}.fam -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/our_samples/ibis/${file}_7
	~/ibis/ibis ${file}.bed ${file}_map.bim ${file}.fam -min_l 5 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/our_samples/ibis/${file}_5
	~/ibis/ibis ${file}.bed ${file}_map.bim ${file}.fam -min_l 0 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/our_samples/ibis/${file}_0
done < populations

# Run IBIS. Which parameters to include?
# min_l = minimum length for acceptable segments to output. I run it once 7cM (default) and once with 5cM (I have also run it with 1cM to check differences).
# mt = minimum number of markers required for a segment to be printed (default is 436). 
# mt2 = minimum number of markers required for an IBD2 segment to be printed (default is 186).
# er = acceptable error rate in a segment before considering it false (default is 0.004).
# printCoef = Have the program print a .coef file and a .incoef (with -hbd) file.
# f = specify output prefix
```

3. Check results of IBIS in R
```
# Executable in R. Based on script IBIS_testing.R in ~/RauteGenetics, created on Wednesday 19 July 2023.
# Script designed to look at files created by IBIS software to examine kinship based on IBD


######################################################################################################################
                                    Import libraries & set working directory
######################################################################################################################

# Import libraries 
library(reshape2)
library(igraph)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(ggraph)

# set working directory

setwd("/Users/inezd/RauteGenetics")

######################################################################################################################
                                    Write functions necessary to adapt files
######################################################################################################################


#### 1. Merge with samples_origins

change_and_merge <- function(x){
  
  names(x) <- c("sample1", "sample2", "kinship_coefficient", "IBD2_fraction", "segment_count", "degree_of_relatedness")
  
  # Change individual names
  x$sample1 <- gsub("^[^:]*:", "", x$sample1)
  x$sample2 <- gsub("^[^:]*:", "", x$sample2)
  
  # return x
  return(x)
}


#### 2. create new selection of individuals who only shared relatedness of PI_HAT < 0.35

create_sub_file <- function(x){
  x <- x %>%
    arrange(sample1, sample2, desc(kinship_coefficient))
  y <- x[,c(1,2,3)]
  names(y)[1:2] <- c("ID1", "ID2")
  y$pair <- paste(pmin(y$ID1, y$ID2), pmax(y$ID1, y$ID2), sep = "-")
  return(y)
}


#### 3. Create file with the right kinship values for all individuals

create_ind_file <- function(x){
  y <- gather(x, ID, type, ID1:ID2, factor_key=TRUE)
  y <- y[,c(1,4)]
  y <- y %>%
    arrange(desc(kinship_coefficient), type)
  return(y)
}


#### 4. Create function toLoop until no kinship coefficient not higher than certain value

loop_kinship <- function(ind_file, coef_file){
  while (any(ind_file$kinship_coefficient > 0.25)) { # Change this value per our input
    # Find the individual (type) with the highest PI_HAT value
    max_kin <- max(ind_file$kinship_coefficient)
    max_individual <- ind_file$type[which(ind_file$kinship_coefficient == max_kin)[1]]
    
    # Delete the individual from df1 and update df2
    coef_file <- coef_file[!(coef_file$ID1 == max_individual | coef_file$ID2 == max_individual), ]
    ind_file <- gather(coef_file, ID, type, ID1:ID2, factor_key=TRUE)
    ind_file <- ind_file[,c(1,4)]
  }
  # who are the unrelated individuals in this case?
  unrelate_ind <- unique(c(coef_file$ID1, coef_file$ID2)) # if kinship coefficient is 0.25, this gives 34 individuals
  unrelate_ind <- unrelate_ind[order(unrelate_ind, decreasing = TRUE)]
  
  # return file
  return (unrelate_ind)
}


#### 5. Check individuals and save output

check_ind <- function(unrelate_ind, coef_file){
  y <- as.data.frame(t(combn(unrelate_ind, 2)))
  names(y)[1:2] <- c("ID1", "ID2")
  y$pair <- paste(pmin(y$ID1, y$ID2), pmax(y$ID1, y$ID2), sep = "-")
  y <- merge(y, coef_file[,3:4], by = 'pair')
  y <- y[order(y$ID1, y$ID2),]
  print(mean(y$kinship_coefficient))
  print(max(check_individuals$kinship_coefficient))
  
  # return output file
  return(y)
  
}


#### 6. Create a datafile ready for saving

save_file <- function(coef_subset, name){
  coef_subset <- data.frame(FID = name, IID = coef_subset)
}


#### 7. Combine all functions into one new function for compactness:

determine_relatedness <- function(coef_file, population){
  coef_file <- change_and_merge(coef_file)
  coef_file_sub <- create_sub_file(coef_file)
  coef_file_ids <- create_ind_file(coef_file_sub)
  coef_file_subset <- loop_kinship(coef_file_ids, coef_file_sub)
  check_individuals <- check_ind(coef_file_subset, coef_file_sub)
  coef_file_subset <- save_file(coef_file_subset, population)
  return(coef_file_subset)
}



######################################################################################################################
                                              Apply functions to Raute and controls data
######################################################################################################################


### Check effect of using different centimorgans

# Import data
Raute7_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/Raute_ibis_7.coef", col_names = T))
Raute5_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/Raute_ibis_5.coef", col_names = T))
Raute0_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/Raute_ibis_0.coef", col_names = T))
Controls_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/controls_ibis_7.coef", col_names = T)) # we use the 7 cM dataset here only


# Apply functions
Raute7_subset <- determine_relatedness(Raute7_coef, "Raute") # 34 individuals  > can just use this one
Raute5_subset <- determine_relatedness(Raute5_coef, "Raute") # 34 individuals
Raute0_subset <- determine_relatedness(Raute0_coef, "Raute") # 32 individuals
Controls_subset <- determine_relatedness(Controls_coef, "Dailekh") # 34 individuals


### Determine smaller amount of individuals to use to check whether the allele frequency analyses improve

# We have to change the loop:
loop_kinship2 <- function(ind_file, coef_file){
  while (any(ind_file$kinship_coefficient > 0.16)) { # we're putting 0.12, because its much lower and will select fewer individuals
    # Find the individual (type) with the highest PI_HAT value
    max_kin <- max(ind_file$kinship_coefficient)
    max_individual <- ind_file$type[which(ind_file$kinship_coefficient == max_kin)[1]]
    
    # Delete the individual from df1 and update df2
    coef_file <- coef_file[!(coef_file$ID1 == max_individual | coef_file$ID2 == max_individual), ]
    ind_file <- gather(coef_file, ID, type, ID1:ID2, factor_key=TRUE)
    ind_file <- ind_file[,c(1,4)]
  }
  # who are the unrelated individuals in this case?
  unrelate_ind <- unique(c(coef_file$ID1, coef_file$ID2)) # if kinship coefficient is 0.25, this gives 34 individuals
  unrelate_ind <- unrelate_ind[order(unrelate_ind, decreasing = TRUE)]
  
  # return file
  return (unrelate_ind)
}

# We now will perform the functions separately:
Raute7_coef <- change_and_merge(Raute7_coef)
Raute7_coef_sub <- create_sub_file(Raute7_coef)
Raute7_coef_ids <- create_ind_file(Raute7_coef_sub)
Raute7_coef_subset <- loop_kinship2(Raute7_coef_ids, Raute7_coef_sub)
check_individuals <- check_ind(Raute7_coef_subset, Raute7_coef_sub)
Raute7_coef_subset <- save_file(Raute7_coef_subset, "Raute")

# Write files: this is for a subset of fewer Raute
write.table(Raute7_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/Raute_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)
write.table(Raute7_coef_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/Raute_subset_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)
write.table(Controls_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/Controls_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)

```

4. Remove the unrelated individuals from the dataset (for now we are ignoring the second individual list we created with a subset of Raute)
```
cd ~/plink/our_samples/ibis/
mkdir ~/plink/our_samples/unrelated_samples

# Keep 34 unrelated Raute
plink --bfile raute_ibis --allow-no-sex \
        --keep ./Unrelated/Raute_unrelatedness.txt \
        --make-bed \
        --out ~/plink/our_samples/unrelated_samples/unrelated_raute

# Keep unrelated controls
plink --bfile controls_ibis --allow-no-sex \
        --keep ./Unrelated/Controls_unrelatedness.txt \
        --make-bed \
        --out ~/plink/our_samples/unrelated_samples/unrelated_controls
```

NOTE: I upload all files from my computer to my server using scp, hence the file names that I save in R and those I use on the VM do not always overlap 


