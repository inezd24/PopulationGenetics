# Genetic structure analyses with flashpca and ADMIXTURE
Here, I show how to use the programs flashpca and ADMIXTURE to examine genetic structure in populations. We use R to visualize the output from flashpca and PONG to visualize the output from ADMIXTURE. 

### Flashpca for pca
Why flashpca? Allows you to look at a lot of different PCs, but also ease of use of program. See https://github.com/gabraham/flashpca for more information. 

1. Install from source
```
git clone git://github.com/gabraham/flashpca

#The Makefile contains three variables that need to be set according to where you have installed the Eigen headers and Boost headers and libraries on your system. The default values for these are:
# EIGEN_INC=/usr/local/include/eigen
# BOOST_INC=/usr/local/include/boost
# BOOST_LIB=/usr/local/lib
# SPECTRA_INC=spectra
# If your system has these libraries and header files in those locations, you can simply run make:
cd flashpca
make all
```

2. Prepare data files for different pcas: first a dataset just with all samples, filtered for maf
```
# For entire dataset (1342 samples)
plink --bfile ~/plink/modified_samples/pheno_ready --allow-no-sex \
	--extract ~/plink/modified_samples/all_pops.prune.in \
	--maf 0.05 \
	--make-bed \
	--out entire_network
```

3. Create network with a smaller subset of Raute (only 5), but with all other samples (1313)
```
# First we make a list of the individuals who are in the populations
sort -n -k 1 ~/plink/modified_samples/pheno_ready.fam > all_individuals.txt

# Keep the five individuals of the Raute and all other individuals
grep -v "Raute" all_individuals.txt > raute5_individuals.txt
grep "8378.CEL" all_individuals.txt  >> raute5.txt
grep "8301R.CEL" all_individuals.txt  >> raute5.txt
grep "8234.CEL" all_individuals.txt  >> raute5.txt
grep "8260R.CEL" all_individuals.txt  >> raute5.txt
grep "8235.CEL" all_individuals.txt  >> raute5.txt
cat raute5.txt >> raute5_individuals.txt
rm raute5.txt

# Turn into plink files
plink --bfile ~/plink/modified_samples/pheno_ready --allow-no-sex \
	--keep raute5_individuals.txt \
        --extract ~/plink/modified_samples/all_pops.prune.in \
        --maf 0.05 \
        --make-bed \
        --out raute_5
```

4. Create network with reduced number of individuals for each population included
```
# Select a maximum of 10 individuals per population
awk -v max_lines_per_value=10 '{
    if (++count[$1] <= max_lines_per_value) {
        print > "reduced_network.txt"
    }
}' all_individuals.txt

# But we want the Raute to be the same five individuals only, just in case
grep -v "Raute" reduced_network.txt > reduced_network2.txt
cat raute5.txt >> reduced_network2.txt
rm reduced_network.txt
mv reduced_network2.txt reduced_network.txt

# Turn into plink files
plink --bfile ~/plink/modified_samples/pheno_ready --allow-no-sex \
        --keep reduced_network.txt \
        --extract ~/plink/modified_samples/all_pops.prune.in \
        --maf 0.05 \
        --make-bed \
        --out reduced_network
```

5. Use pruned data to run flashpca
```
# Make a directory with pca output
mkdir pca_output
cd ~/plink/pca_output

# Perform flashpca on each set of data you want to examine (with 20 pcs: -d 20)
~/flashpca_x86-64 --bfile ~/plink/modified_samples/entire_network -d 20 -f _entirenet.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/reduced_network -d 20 -f _reducednet.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/raute_5 -d 20 -f _raute5.txt
```

### Examine and visualise in R

```
######################################################################################################################
############################################### start of script ######################################################

# Script is executable in R. Called flashpca.R in ~/RauteGenetics.
# Last updated on 1 November 2023. 
# Developed to visualise flashpca output in pca plots

######################################################################################################################
                                            Import libraries & files
######################################################################################################################

# Import libraries
library(readr)
library(readxl)
library(ggplot2)
library(purrr)
library(stringr)
library(ggnewscale)
library(tidyverse)
library(colorblindr)
library(plyr)

# Import files
pve_entirenet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_entirenet.txt", col_names = F)
pcs_entirenet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_entirenet.txt")
pve_reducednet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_reducednet.txt", col_names = F)
pcs_reducednet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reducednet.txt")
pve_raute5 <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_raute5.txt", col_names = F)
pcs_raute5 <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_raute5.txt")
population_data <- read_excel("/Users/inezd/PHD/Data/Population_data_information.xlsx") # This is a file with information about each population

######################################################################################################################
                                                Write functions
######################################################################################################################


# 1. Function to change evec files
change_pcs <- function(pcs_file, populations_information){
  
  # First add the IID to the file
  populations_information$IID <- populations_information$Individual
  
  # Merge with information about the population to which the individual belongs
  pcs_file <- merge(pcs_file, populations_information, by = "IID", all.x=T) #populations information is already structure
  
  # Order dataframe according to continent and then population (the combined one)
  pcs_file <- pcs_file[order(pcs_file$Region, pcs_file$Population),]
  
  # Create temporary file to list numbers for the different shapes
  shape_palette <- pcs_file %>% dplyr::group_by(Region, Population) %>% dplyr::summarise(n=n())
  shape_palette$shape <- ave(shape_palette$n,                 # Create numbering variable
                             shape_palette$Region,
                             FUN = seq_along)
  
  # Order again
  shape_palette <- shape_palette[order(shape_palette$Region, shape_palette$Population),]
  
  # Merge with pcs_file
  pcs_file <- merge(pcs_file, shape_palette, by = c("Region", 'Population'), all.x=T)
  
}

# 2. Function to make fit for splitting
split_pcs <- function(x){
  x$shape <- as.factor(as.character(x$shape))
  x$Population <- as.character(x$Population)
  y <- split(x, x$Region)
  return(y)
}

# 3. Create function o plot first two pcs
plot_PC12 <- function(split_df, color_list, shape_list, col_nr, size_points, PC_A, PC_B,
                      perc_PCA, perc_PCB){
  
  #Define the variables
  colors = color_list
  shapes = shape_list
  ncol = col_nr
  size = size_points
  lab_x = paste('PC', PC_A, ' (', perc_PCA, '%)')
  lab_y = paste('PC', PC_B, ' (', perc_PCB, '%)')
  
  # use ggplot to plot
  ggplot(mapping = aes(x = PC1, y = PC2)) +
    purrr::imap(split_df, function(x, y) {
      
      # Get Labels
      labels <- x[c("shape", "Population")] %>% 
        distinct(shape, Population) %>% 
        deframe()
      
      # Get order
      list(
        geom_point(data = x, aes(shape = shape), color = colors[[y]], size = 3),
        stat_ellipse(data = x, color = colors[[y]]),
        scale_shape_manual(values = shapes, labels = labels, name = y, guide = guide_legend(ncol = ncol, override.aes = list(size = size))),
        theme_light(),
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.title.x.top = element_blank(), 
              axis.text.x.top = element_blank(), 
              axis.title.y.right = element_blank(), 
              axis.text.y.right = element_blank(),
              text = element_text(size = 20),
              panel.background = element_rect(fill='transparent'), 
              plot.background = element_rect(fill='transparent', color=NA), 
              legend.background = element_rect(fill='transparent')),
        labs(x =lab_x, 
             y = lab_y),
        new_scale("shape")
      )
    })
}


# 4. Create function to plot first third and fourth pcs
plot_PC34 <- function(split_df, color_list, shape_list, col_nr, size_points, PC_A, PC_B,
                      perc_PCA, perc_PCB){
  
  #Define the variables
  colors = color_list
  shapes = shape_list
  ncol = col_nr
  size = size_points
  lab_x = paste('PC', PC_A, ' (', perc_PCA, '%)')
  lab_y = paste('PC', PC_B, ' (', perc_PCB, '%)')
  
  # use ggplot to plot
  ggplot(mapping = aes(x = PC3, y = PC4)) +
    purrr::imap(split_df, function(x, y) {
      
      # Get Labels
      labels <- x[c("shape", "Population")] %>% 
        distinct(shape, Population) %>% 
        deframe()
      
      # Get order
      list(
        geom_point(data = x, aes(shape = shape), color = colors[[y]], size = 3),
        stat_ellipse(data = x, color = colors[[y]]),
        scale_shape_manual(values = shapes, labels = labels, name = y, guide = guide_legend(ncol = ncol, override.aes = list(size = size))),
        theme_light(),
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.title.x.top = element_blank(), 
              axis.text.x.top = element_blank(), 
              axis.title.y.right = element_blank(), 
              axis.text.y.right = element_blank(),
              text = element_text(size = 20),
              panel.background = element_rect(fill='transparent'), 
              plot.background = element_rect(fill='transparent', color=NA), 
              legend.background = element_rect(fill='transparent')),
        labs(x =lab_x, 
             y = lab_y),
        new_scale("shape")
      )
    })
}

######################################################################################################################
                                                Apply functions
######################################################################################################################

# Change files
pcs_entirenet <- change_pcs(pcs_entirenet, populations_list)
pcs_entirenet <- unique(pcs_entirenet)
pcs_reducednet <- change_pcs(pcs_reducednet, populations_list)
pcs_reducednet <- unique(pcs_reducednet)
pcs_raute5 <- change_pcs(pcs_raute5, populations_list)
pcs_raute5 <- unique(pcs_raute5)

# Split
pcs_entire_split <- split_pcs(pcs_entirenet)
pcs_reduced_split <- split_pcs(pcs_reducednet)
pcs_raute5_split <- split_pcs(pcs_raute5)

# Change eigenvalues
pve_entirenet$X1 <- round(pve_entirenet$X1*100,2)
names(pve_entirenet) <- 'Eigenvalue'
pve_reducednet$X1 <- round(pve_reducednet$X1*100,2)
names(pve_reducednet) <- 'Eigenvalue'
pve_raute5$X1 <- round(pve_raute5$X1*100,2)
names(pve_raute5) <- 'Eigenvalue'

######################################################################################################################
                                                    Visualize
######################################################################################################################

# Set shapes for all samples
shapes <- setNames(c(seq(1:20)),
                   c(seq(1:20)))
ncol <- 3
size = 3

# Entire network
colors <- setNames(c("#ee8866", "#77aadd", "#eedd88", "#ffaabb","#99ddff", "#838B83",  
                     "#836FFF", "#8B6914", "#ffae34", "#762A83", "#000000"), names(pcs_entire_split))
pcs_entirenet_12 <- plot_PC12(pcs_entire_split, colors, shapes, 3, 3, 1, 2, 5.04, 2.18)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_entirenet_12.png",pcs_entirenet_12, bg='white', height = 18, width = 25)
pcs_entirenet_34 <- plot_PC12(pcs_entire_split, colors, shapes, 3, 3, 1, 2, 1.21, 1.20)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_entirenet_34.png",pcs_entirenet_34, bg='transparent', height = 18, width = 25)

# Reduced network
colors <- setNames(c("#ee8866", "#77aadd", "#eedd88", "#ffaabb","#99ddff", "#838B83",  
                     "#836FFF", "#8B6914", "#ffae34", "#762A83", "#000000"), names(pcs_reduced_split))
pcs_reducednet_12 <- plot_PC12(pcs_reduced_split, colors, shapes, 3, 3, 1, 2, 4.30, 2.07)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reducednet_12.png",pcs_reducednet_12, bg='transparent', height = 18, width = 25)
pcs_reducednet_34 <- plot_PC34(pcs_reduced_split, colors, shapes, 3, 3, 1, 2, 1.56, 0.96)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reducednet_34.png",pcs_reducednet_34, bg='transparent', height = 18, width = 25)

# Entire network, but only 5 Raute
colors <- setNames(c("#ee8866", "#77aadd", "#eedd88", "#ffaabb","#99ddff", "#838B83",  
                     "#836FFF", "#8B6914", "#ffae34", "#762A83", "#000000"), names(pcs_raute5_split))
pcs_raute5_12 <- plot_PC12(pcs_raute5_split, colors, shapes, 3, 3, 1, 2, 4.30, 2.07)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reducednet_12.png",pcs_raute5_12, bg='transparent', height = 18, width = 25)
pcs_raute5_34 <- plot_PC34(pcs_raute5_split, colors, shapes, 3, 3, 1, 2, 1.56, 0.96)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reducednet_34.png",pcs_raute5_34, bg='transparent', height = 18, width = 25)

################################################ end of script #######################################################
######################################################################################################################
```

### ADMIXTURE and PONG

1. Download ADMIXTURE
```
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xvzf admixture_linux-1.3.0.tar.gz
```

2. Run ADMIXTURE (note: this still needs to be redone and is currently run on an old file). Currently, this lets ADMIXTURE run once for Ks 3 to 12. We need to do at least three different runs and give them different seeds. 
```
#create necessary files
BED_FILE=~/plink/redo/pcas/maf_05.bed #this file no longer should be used or exists in that location
OUTDIR=~/plink/clean_admixture
mkdir -p $OUTDIR

#loop through multiple k-values
for K in {3..12}; do
	CMD=cd $OUTDIR; admixture --cv $BED_FILE $K
	echo $CMD
	# sbatch -c 8 --mem 12000 --wrap="$CMD"
done
```

3. Make the connection between the server and your laptop (you should open a new shell and run it from your laptop, not within the server).
```
ssh -L 127.0.0.1:4000:localhost:4000 ubuntu@123.45.67.89 # The latter part should have the IP address of the server instead of the 1-9 sequence.
```

4. Install PONG
```
sudo pip install pong
```

5. Prepare files in PONG format. See https://github.com/ramachandran-lab/pong/tree/master/example_data_1kG-p3 for how this data should look. You should have:
```
ind2pop.txt
pong_filemap
pop_order.txt
pop_order_expandednames.txt # optional
```

7. Run PONG
```
# Path to files
PONG=~/plink/clean_admixture/pong_filemap
IND=~/plink/clean_admixture/admixture_indlist2.txt
ORDER=~/plink/clean_admixture/admixture_orderlist.txt

# execute
pong -m $PONG -i $IND -n  $ORDER
```




