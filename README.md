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

2. Use pruned data to run flashpca
```
# Working directory
mkdir ~/plink/pca_output
cd ~/plink/pca_output

# Perform flashpca
~/flashpca_x86-64 --bfile ~/plink/modified_samples/all_pops_LD -d 20 -f _allpopsLD.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/reduced_network -d 20 -f _reduced.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/raute_subset -d 20 -f _raute_subset.txt
#~/flashpca_x86-64 --bfile ~/plink/modified_samples/regional_network -d 20 -f _regionalnet.txt
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
pve_entirenet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_allpopsLD.txt", col_names = F)
pcs_entirenet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_allpopsLD.txt")
pve_reducednet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_reduced.txt", col_names = F)
pcs_reducednet <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_reduced.txt")
pve_raute_subset <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pve_raute_subset.txt", col_names = F)
pcs_raute_subset <- read_table("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_raute_subset.txt")
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

# I also have the same functions for checking other PC combinations, but you can easily adapt these by change the x and y in the ggplot mapping.


######################################################################################################################
                                                Apply functions
######################################################################################################################

# Change files
pcs_entirenet <- change_pcs(pcs_entirenet, populations_list)
pcs_entirenet <- unique(pcs_entirenet)
pcs_reducednet <- change_pcs(pcs_reducednet, populations_list)
pcs_reducednet <- unique(pcs_reducednet)
pcs_raute_subset <- change_pcs(pcs_raute_subset, populations_list)
pcs_raute_subset <- unique(pcs_raute_subset)

# Split
pcs_entire_split <- split_pcs(pcs_entirenet)
pcs_reduced_split <- split_pcs(pcs_reducednet)
pcs_raute_subset_split <- split_pcs(pcs_raute_subset)

# Change eigenvalues
pve_entirenet$X1 <- round(pve_entirenet$X1*100,2)
names(pve_entirenet) <- 'Eigenvalue'
pve_reducednet$X1 <- round(pve_reducednet$X1*100,2)
names(pve_reducednet) <- 'Eigenvalue'
pve_raute_subset$X1 <- round(pve_raute_subset$X1*100,2)
names(pve_raute_subset) <- 'Eigenvalue'

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

# Entire network, but only few Raute
colors <- setNames(c("#ee8866", "#77aadd", "#eedd88", "#ffaabb","#99ddff", "#838B83",  
                     "#836FFF", "#8B6914", "#ffae34", "#762A83", "#000000"), names(pcs_raute5_split))
pcs_raute_subset_12 <- plot_PC12(pcs_raute_subset_split, colors, shapes, 3, 3, 1, 2, 4.30, 2.07)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_raute_subset_12.png",pcs_raute_subset_12, bg='transparent', height = 18, width = 25)
pcs_raute_subset_34 <- plot_PC34(pcs_raute_subset_split, colors, shapes, 3, 3, 1, 2, 1.56, 0.96)
ggsave("/Users/inezd/PHD/Genetic_Data/PCA_output/pcs_raute_subset_34.png",pcs_raute_subset_34, bg='transparent', height = 18, width = 25)

################################################ end of script #######################################################
######################################################################################################################
```

### Fst

We can compute Fst with different programs: PLINK, EIGENSOFT, and also in R via Stamppp. This requires a rewritten raw file. Here we explore all options. 

1. Make and move to working directory
```
mkdir ~/plink/Fst/
cd ~/plink/Fst/
```

2. Compute with eigensoft
```
# Remake the .fam file where the population name is the phenotype (for the conversion to EIGENSTRAT
awk '{print $1,$2,$3,$4,$5,1}' ~/plink/modified_samples/raute_subset.fam > raute_subset_altered.fam

# Make a file that has all individuals in it
awk '{print $2}' ~/plink/modified_samples/raute_subset.fam | sort | uniq > raute_subset_individuals.txt

# Create first .par file
rm par.file.eig
echo genotypename: /home/ubuntu/plink/modified_samples/raute_subset.bed > par.file.eig
echo snpname: /home/ubuntu/plink/modified_samples/raute_subset.bim >> par.file.eig
echo indivname: raute_subset_altered.fam >> par.file.eig
echo outputformat: EIGENSTRAT >> par.file.eig
echo genotypeoutname: raute_subset_fst.geno >> par.file.eig
echo snpoutname: raute_subset_fst.snp >> par.file.eig
echo indivoutname: raute_subset_fst.ind >> par.file.eig
echo poplistname: raute_subset_individuals.txt >> par.file.eig

echo "File par.file.eig created"

# Create second .par file
rm par.file.fst
echo genotypename: ./raute_subset_fst.geno > par.file.fst
echo snpname: ./raute_subset_fst.snp >> par.file.fst
echo indivname: ./raute_subset_fst.ind >> par.file.fst
echo fstonly: YES >> par.file.fst
echo inbreed: YES >> par.file.fst
echo fstdetailsname: ./raute_subset_fst_info.out >> par.file.fst
echo phylipoutname: ./raute_subset_phyl_fst >> par.file.fst
echo evecoutname: ./raute_subset_fst.evec >> par.file.fst
echo evaloutname: ./raute_subset_fst.eval >> par.file.fst

echo "File par.file.fst created"

# Convert data
~/eigensoft/EIG-7.2.1/bin/convertf -p ./par.file.eig > ./out_convert.log
echo "Data converted"

# Calculate Fst
~/eigensoft/EIG-7.2.1/bin/smartpca -p ./par.file.fst >./out_fst.log
echo "Fst calculated"
```


3. Compute with PLINK
```
echo "Start with PLINK analysis"

### For Raute only
echo "First Raute only"

# If for Raute only: Create file for raute
grep "Raute" ~/plink/modified_samples/all_pops_LD.fam > raute_clusters.txt
awk '{print $1, $2}' raute_clusters.txt > raute_fst.txt
rm raute_clusters.txt
awk '{$3 = $2} 1' raute_fst.txt > raute_clusters.txt
rm raute_fst.txt

# Compute Fst with pruned and filtered data
plink2 --bfile ~/plink/modified_samples/all_pops_LD --allow-no-sex \
        --within raute_clusters.txt \
	--fst CATPHENO method=hudson blocksize=100 \
        --out raute34_fst

### For all samples
echo "Now for all samples"

awk '{print $1, $2}' ~/plink/modified_samples/raute_subset.fam > raute_subset_pops.txt
awk '{$3 = $1} 1' raute_subset_pops.txt > raute_subset_popclusters.txt
rm raute_subset_pops.txt

# Perform PLINK Fst
plink2 --bfile ~/plink/modified_samples/raute_subset --allow-no-sex \
        --within raute_subset_popclusters.txt --fst CATPHENO method=wc \
        --out raute_subset_fst
```


3. Lastly recode to raw files for STaMPP analysis in R
```
# First with all populations
plink --bfile ~/plink/modified_samples/all_pops_LD --allow-no-sex \
	--recode A \
	--out Fst_stampp

echo "Done for all populations"

# Then only with Raute
plink --bfile ~/plink/modified_samples/raute_subset --allow-no-sex \
        --recode A \
        --out Fst_stampp_subset

echo "Done for only Raute"
```

4. Perform visualisation in R
```
### Visualise PLINK Fst

fst_pops <- readr::read_table("/Users/inezd/PHD/Genetic_Data/F_statistics/Fst/raute_subset_fst.fst.summary", 
                                col_names=c("pop1", "pop2", "fst"),
                                col_types='ccc')
fst_pops <- fst_pops[-1,]
fst_pops$fst <- as.numeric(fst_pops$fst)
hist(fst_pops$fst) # Differences are minuscule
fst_all_plot <- ggplot(fst_pops,aes(x=pop1, y=pop2, fill=fst))+
  geom_point(shape=21, size = 4)+
  scale_fill_gradient( low="cyan3",high="darkorchid") +
  theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 13),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.background = element_rect(fill='transparent'))
fst_all_plot


### Visualise sTaMPP Fst

# Import files
PopGentoLight=read.PLINK(file = "/Users/inezd/PHD/Genetic_Data/F_statistics/Fst_stampp.raw",
                         parallel = F)

# Convert to gen light object
popgenlight=stamppConvert(PopGentoLight, "genlight")

# Calculate pairwise Fst
pair_fst=stamppFst(popgenlight, 100, 95, 4)
melt_fst=melt(pair_fst)

# Plot the phylogenies
tree=nj(as.dist(pair_fst))

plot(tree)
plot.phylo(tree, "phylogram", edge.width=2, font=3, tip.color='black', cex=1.7)

# Observe in a heatmap
heatmap_fst <- ggplot(data = melt_fst, 
                      aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "red", name="FST")  + 
  ggtitle(expression(atop("Pairwise FST, WC (1984)", atop(italic("N = 143, L = 3,759"), ""))))+
  labs( x = "Sampling Site", y = "Sampling Site") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + 
  coord_fixed()
ggsave("./PCA_Plots/Redo/heatmap_fst.png",heatmap_fst, bg='transparent', height = 20, width = 20)


### Fst calculation with eigensoft
fst_eigensoft <- readr::read_table("/Users/inezd/Documents/PHD/Data/Oragene/Output/raute_subset_fst_info.out")
names(fst_eigensoft)[1:6] <- c("pop1", "pop2","snpname", "N", "D", 'Ratio')
fst_eigensoft <- fst_eigensoft[,1:6]
fst_eigensoft$Ratio <- as.numeric(fst_eigensoft$Ratio)
fst_eig_mean <- fst_eigensoft %>% 
  group_by(pop1, pop2) %>%
  summarise(N = mean(N), D = mean(D), Ratio = mean(Ratio, na.rm=T))
fst_eig_mean$pop_pair <- paste(pmin(fst_eig_mean$pop1, fst_eig_mean$pop2), 
                               pmax(fst_eig_mean$pop1, fst_eig_mean$pop2), sep = '-')


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
BED_FILE= ~/plink/modified_samples/all_pops_LD.bed # this file is filtered for MAF and pruned for LD
OUTDIR=~/plink/clean_admixture
mkdir -p $OUTDIR

#loop through multiple k-values
for run in {1..3}; do
	for K in {2..18}; do
		seed=$RANDOM
		CMD=cd $OUTDIR; admixture -s $seed --cv $BED_FILE $K  > raute_5.${run}.${K}.log
		echo $CMD

	# Change name
	mv ~/plink/raute_5.${K}.Q ~/plink/raute_5.${K}.${run}.Q
	mv ~/plink/raute_5.${K}.P ~/plink/raute_5.${K}.${run}.P
done
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




