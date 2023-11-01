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

2. Use PLINK 1.9 to filter VCF file and convert to plink format
```
VCF=VCF/GenotypingRaute2.vcf # VCF file with all data
INCL=~/plink/subset_text_files/samples_include.txt # File with the samples that we want to include out of all samples in the VCF (some are doubles for instance)
SEX=~/plink/subset_text_files/update_sex_ids.txt # Sex IDs for some which are ambiguous
plink --vcf $VCF --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-missing-var-ids @:# \
        --keep $INCL \ # only keep the individuals in the list provided
        --update-sex $SEX \ # update the sex ids (1 or 2) with the list provided
        --biallelic-only \ # only keep biallelic SNPs, filters out triallelic ones. 
        --geno 0.05 \ # filters out all variants with missing call rates exceeding the provided value
        --mind 0.1 \ # does the same for samples.
        --hwe 0.000001 midp \ # filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold.
        --make-bed \ # make bed/bim/fam files
        --out ./our_samples/our_samples # name of output file
```

3. Use plink 2 to set variant ids as 'chr_position' (plink 1.9 doesn't do this) and remove duplicates
```
plink2 --bfile ./our_samples/our_samples --double-id --allow-no-sex --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \ 
        --make-bed \ # make bed/bim/fam files
        --rm-dup exclude-all \ # remove duplicates
        --exclude ~/plink/remove_snps_oursamples \ # exclude the snps that you want removed
        --out ./our_samples/our_samples2 # name of output file
```

4. Filter for minor allele frequency and prune for linkage disequilibrium. Pruning creates a subset of markers that are in approximate linkage equilibrium with each other.
```
plink --bfile our_samples2 --double-id --allow-no-sex \
        --maf 0.05 \ # set minor allele frequency to 5%
        --indep-pairwise 200 25 0.5 \ # window size = 200, step size = 25, and r^2 threshold is 0.5. 
        --make-bed \ # make bed/bim/fam files
        --out ./our_samples/our_samples_maf # name of output file

plink --bfile our_samples_maf --double-id --allow-no-sex \ # take the samples also filtered for maf
        --extract ~/plink/our_samples_maf.prune.in \ # extract the pruned data from the previous step to only include those markes in linkage equilibrium
        --make-bed \ # again, make new bed/bim/fam files 
        --out ./our_samples/our_samples_LD # name of output file
# These two data types, 'our_samples_maf' and 'our_samples_LD', are the ones we will continue to learn 
```

### Using IBIS to extract kinship coefficient and determine unrelated individuals in R
See https://github.com/williamslab/ibis/blob/master/README.md for more instructions on how to use IBIS as well as an overview of all the possible parameters. IBIS is aimed at unphased datasets. We use IBIS on the datafiles that are not yet filtered for MAF or pruned for LD. Use of script 'raute_IBIS.sh' from PLINK directory.

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
# Locate files to use
BED=~/plink/our_samples2.bed
FAM=~/plink/our_samples2.fam
BIM=~/plink/our_samples2.bim
SAMP=~/plink/our_samples2

# Insert map into bim.file and locate file
~/ibis/add-map-plink.pl $BIM ~/ibis/genetic_map_GRCh37_chr{1..22}.txt > ~/plink/our_samples2_map.bim
BIM2=~/plink/our_samples2_map.bim

# Run IBIS. Which parameters to include?
# min_l = minimum length for acceptable segments to output. I run it once 7cM (default) and once with 5cM (I have also run it with 1cM to check differences).
# mt = minimum number of markers required for a segment to be printed (default is 436). 
# mt2 = minimum number of markers required for an IBD2 segment to be printed (default is 186).
# er = acceptable error rate in a segment before considering it false (default is 0.004).
# printCoef = Have the program print a .coef file and a .incoef (with -hbd) file.
# f = specify output prefix
~/ibis/ibis $BED $BIM2 $FAM -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/roh/test1 
~/ibis/ibis $BED $BIM2 $FAM -min_l 5 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/roh/test5

```

3. Check results of IBIS in R
```
# Executable in R. Based on script IBIS_testing.R in ~/RauteGenetics, created on Wednesday 19 July 2023.
# Script designed to look at files created by IBIS software to examine kinship based on IBD
# Set working directory
setwd("/Users/inezd/RauteGenetics")

######################################################################################################################
                                    Import libraries & files
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

# Import files and give correct variable names
test7 <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/Output/test1.seg", col_names = F))
names(test7) <- c("sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos", 
                  "IBD_type", "genetic_start_pos", "genetic_end_pos", "genetic_seg_length",
                  "marker_count", "error_count", "error_density")
test7_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/Output/test1.coef", col_names = T))
names(test7_coef) <- c("sample1", "sample2", "kinship_coefficient", "IBD2_fraction", "segment_count", "degree_of_relatedness")


######################################################################################################################
                                    Write functions necessary to adapt files
######################################################################################################################


# 1. Function to adapt files
change_and_merge <- function(x, y){
  
  # Change individual names
  x$sample1 <- gsub("^[^:]*:", "", x$sample1)
  x$sample2 <- gsub("^[^:]*:", "", x$sample2)
  
  # Make temporary file so you don't have to change anything
  tmp <- y
  
  # merge by sample 1
  names(tmp)[2] <- 'sample1'
  x <- merge(x, tmp[,c(2,4)], by = 'sample1', all.x=T)
  names(x)[7] <- 'origin_sample1'
  
  # merge by sample 2
  names(tmp)[2] <- 'sample2'
  x <- merge(x, tmp[,c(2,4)], by = 'sample2', all.x=T)
  names(x)[8] <- 'origin_sample2'
  
  # re-order columns and rows
  x <- x[order(x$sample1, x$sample2),c(2,7,1,8,3:6)]
  rownames(x) <- seq_along(1:nrow(x))
  
  # return x
  return(x)
}

# 2. Function to create new selection of individuals who only shared relatedness of PI_HAT < 0.35. PI_HAT is the kinship coefficient from IBIS based on IBDs
create_sub_file <- function(x){
  x <- x %>%
    arrange(sample1, sample2, desc(kinship_coefficient))
  y <- x[,c(1,3,5)]
  names(y)[1:2] <- c("ID1", "ID2")
  y$pair <- paste(pmin(y$ID1, y$ID2), pmax(y$ID1, y$ID2), sep = "-")
  return(y)
}

# 3. Function to create file with the right kinship values for all individuals
create_ind_file <- function(x){
  y <- gather(x, ID, type, ID1:ID2, factor_key=TRUE)
  y <- y[,c(1,4)]
  y <- y %>%
    arrange(desc(kinship_coefficient), type)
  return(y)
}

# 4a. Function to create function toLoop until no kinship coefficient not higher than certain value
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

# 4b. Alternative option to 'loop_kinship', but doesn't work well on my computer (crashes with certain sample sizes)
network_kinship <- function(coef_file){
  
  # create subset file
  subset <- subset(coef_file, coef_file$kinship_coefficient <= 0.25)
  subset <- subset[order(subset$ID1, subset$ID2),]
  
  # create network
  g <- igraph::graph_from_data_frame(subset[,1:3], directed = FALSE)
  g <- set_edge_attr(g, "weight", value= subset$kinship_coefficient)
  
  # find largest weighted cliques
  groups <- maximal.cliques(g) # this is where it normally crashes at some point for me
  
  # look for mean edge weights per clique
  mean_edge_weights <- sapply(groups, function(clique) {
    subgraph <- induced.subgraph(net, clique)
    mean(E(subgraph)$weight)
    })
  
  # Find the index of the clique with the highest mean edge weight
  min_mean_weight_index <- which.min(mean_edge_weights)
  
  # Get the clique with the highest mean edge weight
  min_mean_weight_clique <- groups[[min_mean_weight_index]]
  
  # Convert to df
  nodes <- V(net)[min_mean_weight_clique]$name
  
  #print the relevant stats
  print(groups)
  print(mean_edge_weights)
  print(min_mean_weight_index)
  print(min_mean_weight_clique)
  return(nodes)
}

# 5. Function to check individuals and save output
check_ind <- function(unrelate_ind, coef_file){
  y <- as.data.frame(t(combn(unrelate_ind, 2)))
  names(y)[1:2] <- c("ID1", "ID2")
  y$pair <- paste(pmin(y$ID1, y$ID2), pmax(y$ID1, y$ID2), sep = "-")
  y <- merge(y, coef_file[,3:4], by = 'pair')
  y <- y[order(y$ID1, y$ID2),]
  print(mean(y$kinship_coefficient))
  
  # return output file
  return(y)
  
}

######################################################################################################################
                                              Apply functions to data
######################################################################################################################

# Apply all functions
test7_coef <- change_and_merge(test7_coef, samples_origins) # step 1
raute7_coef <- subset(test7_coef, origin_sample1 == 'Raute' & origin_sample2 == 'Raute') # only keep cases where both individuals are Raute
raute7_coef_sub <- create_sub_file(raute7_coef) # step 2
raute7_ids <- create_ind_file(raute7_coef_sub) # step 3
raute7_subset <- loop_kinship(raute7_ids, raute7_coef_sub) # step 4a
# or use alternative: raute7_subset2 <- network_kinship(raute7_coef_sub) step 4b
check_individuals <- check_ind(raute7_subset, raute7_coef_sub) # step 5
# or use alternative: check_individuals <- check_ind(raute7_subset2, raute7_coef_sub) # step 5

# Look at distribution of overall kinship 
hist_coef <- ggplot(raute7_coef, aes(x=kinship_coefficient)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  theme_minimal()+
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x = 'Kinship coefficient', y = 'density') +
  theme(text = element_text(size = 20))
ggsave("./stats_plots/hist_coef.png", hist_coef, bg='transparent', height = 17, width = 17)

# Plot distribution of kinship coefficient of included samples (only those 34 individuals)
hist_kinshpcoef_incl_7 <- ggplot(check_individuals, aes(x = kinship_coefficient)) +
  geom_histogram(col = 'dark grey') +
  theme_minimal() +
  labs(x = 'Kinship coefficient', y = 'Frequency')
mean(check_individuals$kinship_coefficient)
median(check_individuals$kinship_coefficient)
max(check_individuals$kinship_coefficient)
min(check_individuals$kinship_coefficient)

# Save file
raute_subset <- data.frame(FID = raute7_subset, IID = raute7_subset)
write.table(raute_subset,'./OutputFiles/raute_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)

# The latter output we use for omitting related individuals in plink

```

4. Run PLINK with unrelated (unrelated here equals the PI_HAT threshold we set in R) individuals.
```
SAMP=~/plink/our_samples/our_samples2
raute_keep=~/plink/subset_text_files/raute_unrelatedness_34.txt
plink --bfile $SAMP --double-id --allow-no-sex \
        --keep $raute_keep \
        --make-bed \
        --out ~/plink/our_samples/unrelated_raute
```

NOTE: I upload all files from my computer to my server using scp, hence the file names that I save in R and those I use on the VM do not always overlap 


