# Merging SNP array data with PLINK
Here, I merge data from different types of arrays, e.g. different Illumina arrays and the Human Origins array, to produce a large dataset. Almost all data I use is public and can be obtained via https://reich.hms.harvard.edu/datasets. Some datasets I use for my analyses are not public and therefore I do not show how to download any datasets, but only call them by first author and year. I also run IBIS on all of these datasets to make sure that the individuals included are minimally related or at least in the same way as the PI_HAT threshold I set for my own samples. 

### Filter all datasets (even if from different arrays) and check for related individuals

1. Create directory for data from other labs.
```
mkdir ~/plink/data_others
# Go into directory and download all data that you want to add using for instance 'wget <link>' if public. If private, use way specified by authors/committee.
```

2. Download EIGENSOFT to use function 'convertf'
```
cd; wget https://github.com/DReichLab/EIG/archive/v7.2.1.tar.gz
gunzip v7.2.1.tar.gz
mv v7.2.1 eigensoft
cd eigensoft/src
sudo make all
sudo make install
```

3. Create a par file for each dataset that needs to be converted in the following format:
```
genotypename: input genotype file
snpname:      input snp file
indivname:    input indiv file 
outputformat:    ANCESTRYMAP,  EIGENSTRAT, PED, PACKEDPED or PACKEDANCESTRYMAP
                 (Default is PACKEDANCESTRYMAP.)
genotypeoutname: output genotype file
snpoutname:      output snp file
indivoutname:    output indiv file
```

4. Use convertf function to convert between data formats if necessary
```
~/eigensoft/EIG-7.2.1/bin/convertf -p par.file
```

5. All datasets need to be filtered in the same way
```
# Refer to files
WANG=~/plink/data_others/wang_data/Wang_2021
LAZA=~/plink/data_others/lazaridis_data/Lazaridis_2014
NAKA=~/plink/data_others/nakatsuka_data/Nakatsuka_2017
SKOG=~/plink/data_others/skoglund_data/Skoglund_2016
CHANG=~/plink/data_others/changmai_data/Changmai_2022
JEONG=~/plink/data_others/jeong_data/Jeong_2019
AGTA=~/plink/data_others/agta_data/Agta_2023
ARCI=~/plink/data_others/Arciero_2018
GNEC=~/plink/data_others/gnecchi_ruscone_data/Gnecchi_2017

# New output directory
new_out=~/plink/data_others

# For all other datasets, use a for loop to apply the same filters
for f in $WANG $CHANG $SKOG $JEONG $NAKA $AGTA $GNEC $ARCI
do

	# Extract the filename without path
	filename=$(basename "$f")

	# Append "_clean" to the filename
	new_filename="clean_${filename##*.}"

	# Construct the new output path
	output_path="${new_out}/${new_filename}"

	# Perform plink
	plink2 --bfile $f --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
        --geno 0.05 \
        --mind 0.1 \
        --hwe 0.000001 \
        --make-bed \
        --out "$output_path"
done # we only filter for MAF later when merging datasets

# One file has weird chromosome names, hence we filter separately.
plink2 --bfile $LAZA --chr-set 90 --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
        --geno 0.05 \
        --mind 0.1 \
        --hwe 0.000001 \
        --make-bed \
        --out ~/plink/data_others/clean_Lazaridis_2014

```

6. We have to check if any of these samples still contain any related individuals according to the threshold we set. Add map to all the bim files in the directory so that we can run IBIS.
```
directory_path=~/plink/data_others
EXT=.bim
for f in "$directory_path"/*"$EXT"
do
	filename=$(basename "$f")
	extension="${filename##*.}"
	filename_no_extension="${filename%.*}"
	new_filename="${filename_no_extension}_map.${extension}"
	dir_path=$(dirname "$f")
	output_path="${dir_path}/${new_filename}"
	~/ibis/add-map-plink.pl $f ~/ibis/genetic_map_GRCh37_chr{1..22}.txt > "$output_path"
done
```

7. Run IBIS on each of these datasets.
```
# List of file names (without extensions)
files=("clean_Wang_2021" "clean_Lazaridis_2014" "clean_Jeong_2019" "clean_Skoglund_2016" "clean_Changmai_2022" "clean_Nakatsuka_2017" "clean_Agta_2023" "clean_Arciero_2018")

# Loop through the file names
for file in "${files[@]}"; do
    # Define the input file names
    bed_file="${file}.bed"
    bim_file="${file}_map.bim"
    fam_file="${file}.fam"

    # Run the IBIS program with the specified parameters
    ~/ibis/ibis "$bed_file" "$bim_file" "$fam_file" -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/out

    # Optionally, you can add more commands or processing here for each file
done
```

8. Perform the same action as I did after selecting only unrelated individuals with own samples if necessary for any datasets. For instance:
```
plink2 --bfile ~/plink/data_others/clean_Agta_2023 \
	--keep ~/plink/roh/agtaHG_unrelatedness.txt \ # this is a list of unrelated individuals in the dataset according to chosen PI_HAT value
	--make-bed \
	--out clean_AgtaHunterGatherer_unrl
```


### Merging files from the same array (Human Origins)
Merging files from the same array is relatively straightforward and requires only to create a list of the dataset files you want to merge. 

1. Create a list with files to merge. It should be in order of '.bed .bim .fam'. For instance:
```
cat >> merge_data_new
# dataset1.bed dataset1.bim dataset1.fam
# dataset2.bed dataset2.bim dataset2.fam
# dataset3.bed dataset3.bim dataset3.fam
# datasetx.bed datasetx.bim datasetx.fam
# ctrl+C
```

2. Sort all SNPS per dataset and store these in a different file
```
cd ~/plink/data_others
# Define files > these should be the .bim files that you will use in your merging
files="dataset1.bim dataset2.bim dataset3.bim datasetx.bim"

# Loop through files
for file in ${files}; do

# Check which file we are on
echo $file

# Create a prefix for each file that is to be saved based on name of input file
file_prefix=$(basename "$file" | cut -d. -f1)

# Create output file based on prefix
output_file="${file_prefix}_sorted_SNPs.txt"

# Create new file with all sorted SNPs
awk '{print $2}' $file | sort > "$output_file"

echo "Values from column $file have been extracted to $output_file"

done
```

3. Extract overlapping SNPs between files based.
```
cd ~/plink/data_others
cat *_sorted_SNPs.txt > combined_snps.txt # save all files into one new file
sort combined_snps.txt -o combined_snps.txt # sort the file
uniq -c combined_snps.txt > frequency_snps.txt # create a frequency table out of file. Those with '9' are present in all (9) datasets. 
awk '$1 == 9 {print $2}' frequency_snps.txt > overlapping_snps.txt # only extract those which are present in all datasets and save into new file
```

4. Now merge human origins data
```
cd ~/plink/data_others
raute_unrl=~/plink/our_samples/unrelated_raute # take the unrelated individuals from our sample
plink --bfile $raute_unrl --allow-no-sex \
        --merge-list ~/plink/merging_files/merge_data_new \ # list created in step 1
        --extract ~/plink/data_others/overlapping_snps.txt \ # list created in step 3
        --make-bed \ # make bed bim fam files
        --out ~/plink/modified_samples/HO_samples # output prefix
```

### Merging files from different arrays
If you use datasets that were genotyped on different arrays, there can be several issues that you might encounter. These steps combat some of the issues I encountered and show how dealt with them. 

1. Extract SNP IDs from both .bim files and sort them numerically.
```
BIM_all=~/plink/modified_samples/HO_samples.bim
BIM_arc=~/plink/data_others/clean_Arciero_2018.bim
awk '{print $2}' $BIM_all | sort > clean_samples_sorted.txt
awk '{print $2}' $BIM_arc | sort > clean_arciero_sorted.txt
```

2. Find and extract the IDs common in both files
```
cat clean_samples_sorted.txt clean_arciero_sorted.txt > combined_arc.txt
sort combined_arc.txt -o combined_arc.txt
uniq -c combined_arc.txt > frequency_arc.txt
awk '$1 == 2 {print $2}' frequency_arc.txt >  intersection_arc.txt
```

3. Which variants are double? > as a result of changing the rsid to chr_pos.
```
plink --bfile ~/plink/data_others/clean_Arciero_2018 --write-snplist --out ./arciero_all_snps
cat arciero_all_snps.snplist | sort | uniq -d > arciero_dup_snps.snplist
plink --bfile ~/plink/data_others/clean_Arciero_2018 --exclude ~/plink/merging_files/arciero_dup_snps.snplist --make-bed --out ~/plink/merging_files/arciero_no_dup
```

4. Now extract the snps that are in common with the other dataset
```
plink --bfile ~/plink/modified_samples/HO_samples --allow-no-sex \
	--bmerge ~/plink/merging_files/arciero_no_dup.bed ~/plink/merging_files/arciero_no_dup.bim ~/plink/merging_files/arciero_no_dup.fam \
	--extract ~/plink/merging_files/intersection_arc.txt \
	--make-bed \
	--out ~/plink/merging_files/arciero_ready
```

5. This gives an error message and create a file that ends with '.missnp'. We need to flip the data to take care of these cases.
```
plink --bfile ~/plink/merging_files/arciero_no_dup \
	--flip ~/plink/merging_files/arciero_ready-merge.missnp \
	--make-bed \
	--out ~/plink/modified_samples/flipped_arciero
```

6. Try again to merge > this should succeed
```
plink --bfile ~/plink/modified_samples/HO_samples --allow-no-sex \
        --bmerge ~/plink/modified_samples/flipped_arciero.bed ~/plink/modified_samples/flipped_arciero.bim ~/plink/modified_samples/flipped_arciero.fam \
	--extract ~/plink/merging_files/intersection_arc.txt \
	--indep-pairwise 200 25 0.5 \
	--make-bed \
	--out ~/plink/modified_samples/all_pops
```


























