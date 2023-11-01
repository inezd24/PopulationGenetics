# Genetic structure analyses with flashpca and ADMIXTURE
Here, I show how to use the programs flashpca and ADMIXTURE to examine genetic structure in populations. 

### flashpca for pca
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
mkdir pca_output
cd ~/plink/pca_output

# Perform flashpca on each set of data you want to examine
~/flashpca_x86-64 --bfile ~/plink/modified_samples/entire_network -d 20 -f _entirenet.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/reduced_network -d 20 -f _reducednet.txt
~/flashpca_x86-64 --bfile ~/plink/modified_samples/raute_5 -d 20 -f _raute5.txt
```


