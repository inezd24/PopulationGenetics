# Finestructure 
Pipeline to use the program finestructure with merged genome-wide SNP array data from differented populations. This requires phasing data, which you can do with several programs. I have used BEAGLE (per recommendation of others + ease of use) and you can find the pipeline to do this under the branch 'IBD_analyses'. This software does not like missing data, hence filter for missing genotypes (and SNPs) as well as minor allele frequencies prior to using your data. 

### Preparing files for chromopainter input (test run)
Necessary scripts for this part are try_out_chrom.sh and chr_EM10.sh as well as the perl conversion scripts vcf2cp.pl and convertrecfile.pl. On my server, a lot of the downloading I have not put in a script but done separately. I've tried to add most of it here as otherwise I won't be able to replicate it, but some may be missing (e.g. the prereqs for the GUI took me quite a while to install). 

1. Download, rename, and install finestructure
```
cd; wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.1.1.zip
unzip fs_4.1.1.zip; rm fs_4.1.1.zip
mv fs_4.1.1 finestructure
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
mv ChromoPainterv2 finestructure/
```

4. Download and unzip references files/genetic maps. 
```
git clone https://github.com/adimitromanolakis/geneticMap-GRCh37
gunzip ~/geneticMap-GRCh37/genetic_map_GRCh37_*
```

5. Make test directory and enter directory. 
```
mkdir chr1; cd ./chr1
```

5. Copy BEAGLE phased files into this working directory
```
for chr in {1,5,8,15}; do
cp ~/plink/phased_data/phased_no_LD/single_run/nonLD_phased${chr}.vcf.gz ./
gunzip nonLD_phased${chr}.vcf.gz
done
```

6. Convert these files to chromopainter format via the safe VCF raute provided by finestructure (vcf2cp.pl). See https://www.biostars.org/p/9519621/ or https://github.com/danjlawson/finestructure4 for extra help and the necessary code. 
```
for chr in {1,5,8,15}; do
perl ~/finestructure/vcf2cp.pl nonLD_phased${chr}.vcf nonLD_phased${chr}
done
```

7. Make recombination file with convertrecfile.pl. Do this for a subset of chromosomes only.
```
for chr in {1,5,8,15}; do
~/finestructure/convertrecfile.pl -M hap ~/finestructure/chr1/nonLD_phased${chr}.phase ~/geneticMap-GRCh37/genetic_map_GRCh37_chr${chr}.txt chr${chr}.recombfile 
done
```

8. Run a test with chromopainterv2 for the selected chromosomes. This will be necessary for when we run 
```
for chr in {1,5,8,15}; do
~/finestructure/ChromoPainterv2 -g nonLD_phased${chr}.phase -r chr${chr}.recombfile -t nonLD_phased${chr}.ids -o chr${chr}_EM10 -a 0 0 -i 10 -in -iM
done
```


