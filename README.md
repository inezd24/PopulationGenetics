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
poplistname:	poplist file (which individuals you want to keep)
```


4. Use convertf function to convert between data formats (if necessary). With using convertif, we also update the population list. 
```
cd ~/plink/data_others/lazaridis_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.lazaridis

cd ~/plink/data_others/wang_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.wang

cd ~/plink/data_others/skoglund_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.skoglund

cd ~/plink/data_others/nakatsuka_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.nakatsuka

cd ~/plink/data_others/jeong_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.jeong

cd ~/plink/data_others/changmai_data

~/eigensoft/EIG-7.2.1/bin/convertf -p par.changmai
```


5. We move all newly generated PLINK format files into one directory
```
mkdir ~/plink/data_others/plink_files

# Move all files into this directory
mv ~/plink/data_others/arciero_data/Arciero_2018* ~/plink/data_others/plink_files
mv ~/plink/data_others/agta_data/AgtaHunterGatherer* ~/plink/data_others/plink_files
mv ~/plink/data_others/lazaridis_data/Lazaridis_2014* ~/plink/data_others/plink_files
mv ~/plink/data_others/wang_data/Wang_2021* ~/plink/data_others/plink_files
mv ~/plink/data_others/nakatsuka_data/Nakatsuka_2017* ~/plink/data_others/plink_files
mv ~/plink/data_others/skoglund_data/Skoglund_2016* ~/plink/data_others/plink_files
mv ~/plink/data_others/changmai_data/Changmai_2022* ~/plink/data_others/plink_files
mv ~/plink/data_others/jeong_data/Jeong_2019* ~/plink/data_others/plink_files
```


6. Refer to new files and update IDs and phenotypes
```
# All population files:
WANG=~/plink/data_others/plink_files/Wang_2021
LAZA=~/plink/data_others/plink_files/Lazaridis_2014
NAKA=~/plink/data_others/plink_files/Nakatsuka_2017
SKOG=~/plink/data_others/plink_files/Skoglund_2016
CHANG=~/plink/data_others/plink_files/Changmai_2022
JEONG=~/plink/data_others/plink_files/Jeong_2019
AGTA=~/plink/data_others/plink_files/AgtaHunterGatherer
ARCI=~/plink/data_others/plink_files/Arciero_2018

# Make new directory for files
mkdir /home/ubuntu/plink/data_others/updated_plink_files

# Go into directory we have with files
cd /home/ubuntu/plink/data_others/plink_files

# Path to the folder containing the update files (phenotype and IDs)
update_folder="/home/ubuntu/plink/data_others/fam_files"
plink_folder="."

# Loop through all PLINK files in the specified folder
for f in $WANG $CHANG $SKOG $JEONG $NAKA $AGTA $ARCI; do
	# Extract the filename without path and extension
	filename=$(basename "$f")

	# Construct the paths for phenotype and ID update files
	pheno_file="${update_folder}/update_${filename}_pheno.txt"
	id_file="${update_folder}/update_${filename}_IDS.txt"

	# Replace with your desired output folder
	new_out="/home/ubuntu/plink/data_others/updated_plink_files"

        # Construct the new output path
        output_path="${new_out}/${filename}"

	# Perform PLINK operations with phenotype and ID updates
        plink --bfile "${plink_folder}/${filename}" --allow-extra-chr --allow-no-sex \
	--pheno "${pheno_file}" \
	--update-ids "${id_file}" \
	--make-bed \
	--out "$output_path"
done


## FOR LAZARIDIS DATASET: this one has a weird chromosome (hence we do the --chr-set 90), so we filter it separately. 
filename=$(basename $LAZA)
pheno_file="${update_folder}/update_${filename}_pheno.txt"
id_file="${update_folder}/update_${filename}_IDS.txt"
new_out="/home/ubuntu/plink/data_others/updated_plink_files"
output_path="${new_out}/${filename}"
plink --bfile "${plink_folder}/${filename}" --chr-set 90 --allow-extra-chr --allow-no-sex \
        --pheno "${pheno_file}" \
        --update-ids "${id_file}" \
        --make-bed \
        --out "$output_path"
```


7. Filter each individual population from the different datasets
```
# Enter new directory with all the plink files in it
mkdir ~/plink/data_others/population_files_HO
cd ~/plink/data_others/population_files_HO

# Now, we want to filter them per population:
WANG=~/plink/data_others/updated_plink_files/Wang_2021
LAZA=~/plink/data_others/updated_plink_files/Lazaridis_2014
NAKA=~/plink/data_others/updated_plink_files/Nakatsuka_2017
SKOG=~/plink/data_others/updated_plink_files/Skoglund_2016
CHANG=~/plink/data_others/updated_plink_files/Changmai_2022
JEONG=~/plink/data_others/updated_plink_files/Jeong_2019
AGTA=~/plink/data_others/updated_plink_files/AgtaHunterGatherer
ARCI=~/plink/data_others/updated_plink_files/Arciero_2018

# Filter through each datafile
for f in $WANG $CHANG $SKOG $JEONG $NAKA $AGTA; do

	# Select only per family (population)
	for fam in $(awk '{print $1}' ${f}.fam | sort | uniq); do

		# First print the population
		echo $fam | plink2 --bfile $f --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
		--keep-fam /dev/stdin \
		--geno 0.1 \
		--mind 0.1 \
		--hwe 0.000001 midp \
		--make-bed \
		--out $fam

		# Print the first and last lines of the files
		cat $fam.bim | head -10
		cat $fam.bim | tail -10
	done
done



# Then for Lazaridis only
for fam in $(awk '{print $1}' ${LAZA}.fam | sort | uniq); do

                # First print the population
                echo $fam | plink2 --bfile $LAZA --chr-set 90 --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
                --keep-fam /dev/stdin \
                --geno 0.1 \
                --mind 0.1 \
                --hwe 0.000001 midp \
                --make-bed \
                --out $fam
done

# Now for Arciero, since this data is from a different array, we put it in a different folder so we can merge it separately later. 
mkdir ~/plink/data_others/population_files_ARC
cd ~/plink/data_others/population_files_ARC

for fam in $(awk '{print $1}' ${ARCI}.fam | sort | uniq); do

                # First print the population
                echo $fam | plink2 --bfile $ARCI --allow-extra-chr --chr {1..22} --set-all-var-ids @_# \
                --keep-fam /dev/stdin \
                --geno 0.1 \
                --mind 0.1 \
                --hwe 0.000001 midp \
                --make-bed \
                --out $fam
done
```


8. Peform a relatedness check with IBIS for each population.
```
# New directory
mkdir ~/plink/data_others/ibis_populations

# Go into directory with files
cd ~/plink/data_others/population_files_HO

# Save all filenames
find ~/plink/data_others/population_files_HO/ -type f -name '*.bed' | awk -F/ '{print $NF}' | cut -d. -f1 | sort -u > populations

# Perform IBIS for each population
while IFS= read -r file; do
	~/ibis/add-map-plink.pl ${file}.bim ~/ibis/genetic_map_GRCh37_chr{1..22}.txt > ~/plink/data_others/population_files_HO/${file}_map.bim
	~/ibis/ibis ${file}.bed ${file}_map.bim ${file}.fam -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/data_others/ibis_populations/${file}
done < populations

# Now remove _map.bim files
rm ~/plink/data_others/population_files_HO/*_map.bim

### Now do the same for the populations from Arciero (2018)

cd ~/plink/data_others/population_files_ARC

# Save all filenames
find ~/plink/data_others/population_files_ARC/ -type f -name '*.bed' | awk -F/ '{print $NF}' | cut -d. -f1 | sort -u > populations

# Perform IBIS for each population
while IFS= read -r file; do
	~/ibis/add-map-plink.pl ${file}.bim ~/ibis/genetic_map_GRCh37_chr{1..22}.txt > ~/plink/data_others/population_files_ARC/${file}_map.bim
	~/ibis/ibis ${file}.bed ${file}_map.bim ${file}.fam -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/plink/data_others/ibis_populations/${file}
done < populations

rm ~/plink/data_others/population_files_ARC/*_map.bim
```

9. Check which populations have related individuals based on IBIS results
```
# Check you're into the right folder
cd ~/plink/data_others/ibis_populations

# Find all unique populations in the folder
find ~/plink/data_others/ibis_populations/ -type f -name '*.coef'  | awk -F/ '{print $NF}' | cut -d. -f1 | sort -u > populations

# Check which populations have individuals with higher relatedness
while IFS= read -r file; do
	awk -F'\t' 'NR > 1 && $3 > 0.25 { print FILENAME; exit }' ${file}.coef
done < populations
```


10. Three populations are listed as having related individuals: Kusunda, Agta, and BaYaka. We remove these with R
```
### THIS SCRIPT IS EXECUTABLE IN R TO REMOVE RELATED INDIVIDUALS.
### You need the R script from the previous github branch (filter individuals with PLINK)

# From IBIS analyses, three populations had kin_coef values higher than 0.25: the Agta, BaYaka, and Karitiana

# Import these datafiles:
Agta_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/Agta.coef", col_names = T))
BaYaka_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/BaYaka.coef", col_names = T))
Kusunda_coef <- as.data.frame(read_table("/Users/inezd/PHD/Genetic_Data/IBIS/Kusunda.coef", col_names = T))

# Search for unrelated individuals
Agta_subset <- determine_relatedness(Agta_coef, "Agta") # 65 individuals
BaYaka_subset <- determine_relatedness(BaYaka_coef, "BaYaka") # 12 individuals
Kusunda_subset <- determine_relatedness(Kusunda_coef, "Kusunda") # 8 individuals

# Save as files
write.table(Agta_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/Agta_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)
write.table(BaYaka_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/BaYaka_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)
write.table(Kusunda_subset,'/Users/inezd/PHD/Genetic_Data/IBIS/Unrelated/Kusunda_unrelatedness.txt', sep=" ",col.names = F, row.names = F, quote=F)
```


11. Remove the related individuals from the datasets for the Kusunda, Agta and BaYaka
```
cd ~/plink/data_others/population_files_HO

# Subset Kusunda
plink --bfile Kusunda --allow-no-sex \
        --keep ~/plink/data_others/ibis_populations/Unrelated/Kusunda_unrelatedness.txt \
        --make-bed \
        --out ~/plink/data_others/population_files_HO/Kusunda_unrl

# Subset BaYaka
plink --bfile BaYaka --allow-no-sex \
        --keep ~/plink/data_others/ibis_populations/Unrelated/BaYaka_unrelatedness.txt \
        --make-bed \
        --out ~/plink/data_others/population_files_HO/BaYaka_unrl

# Subset Agta
plink --bfile Agta --allow-no-sex \
        --keep ~/plink/data_others/ibis_populations/Unrelated/Agta_unrelatedness.txt \
        --make-bed \
        --out ~/plink/data_others/population_files_HO/Agta_unrl

# Remove the files with the related individuals
rm ~/plink/data_others/population_files_HO/Agta.*
rm ~/plink/data_others/population_files_HO/BaYaka.*
rm ~/plink/data_others/population_files_HO/Kusunda.*
```


12. Create sorted SNP lists for each population (including the unrelated Raute and our Dailekh control samples)
```
# Working directory
cd ~/plink/data_others/population_files_HO/

# Create new directory for SNPS
mkdir ~/plink/data_others/sorted_SNPs

# Copy Raute and controls to population_files_HO
cp ~/plink/our_samples/unrelated_samples/unrelated_raute* ./
cp ~/plink/our_samples/unrelated_samples/unrelated_controls* ./

# Define files
find ~/plink/data_others/population_files_HO/ -type f -name '*.bim' | awk -F/ '{print $NF}' | sort -u > files

# Loop through files
while IFS= read -r file; do

# Check which file we are on
echo $file

# Create a prefix for each file that is to be saved based on name of input file
file_prefix=$(basename "$file" | cut -d. -f1)
echo $file_prefix

# Create output file based on prefix
output_file="${file_prefix}_sorted_SNPs.txt"
echo $output_file

# Create new file with all sorted SNPs
awk '{print $2}' $file | sort > "$output_file"

echo "Values from column $file have been extracted to $output_file"

# Print number of SNPs
wc -l $output_file

done < files

# Move all files into new directory
mv ./*_sorted_SNPs.txt ~/plink/data_others/sorted_SNPs

# Lastly, we want to remove the Raute and control files from the working directory again.
rm ./unrelated_raute*
rm ./unrelated_controls*
```

13. Now we want to find the overlapping SNPs between all of these populations
```
# Working directory
cd ~/plink/data_others/sorted_SNPs

# Combine all files into one file
cat ./*_sorted_SNPs.txt > combined_snps.txt

# Sort this file
sort combined_snps.txt -o combined_snps.txt

# Find all unique files and sort them by frequency
uniq -c combined_snps.txt > frequency_snps.txt

# Find the number of populations so that we can select only the SNPs which occur in all pops
awk '{if ($1 > max) max = $1} END {print max}' frequency_snps.txt

# This should be 72 (as of Tue 19 December 2023). If the number of populations changes, then this number should be different too in the next cmd. 

# Select only the SNPs which occur in all pops
awk '$1 == 72 {print $2}' frequency_snps.txt > overlapping_snps.txt

# This gives: 158011 (Tue 19 Dec 2023) SNPs which they occur in EACH single population!
wc -l ~/plink/data_others/sorted_SNPs/overlapping_snps.txt
```


14. For each population, select only these 158011 SNPs. This is faster and easier than extracting the SNPs during the merge (I did this previously, and some of the SNPs that will NOT be extracted are duplicates and that causes issues during the merge)
```
# Working directory
cd ~/plink/data_others/population_files_HO/

# Do for all files
find ~/plink/data_others/population_files_HO/ -type f -name '*.bed' | awk -F/ '{print $NF}' | cut -d. -f1 | sort -u > datasets

# Make new directory
mkdir ~/plink/data_others/extracted_SNPs

# Perform IBIS for each population
while IFS= read -r file; do

plink2 --bfile ./$file \
		--extract ~/plink/data_others/sorted_SNPs/overlapping_snps.txt \
		--make-bed \
		--out ~/plink/data_others/extracted_SNPs/$file
done < datasets

# Remove files in the ~/plink/data_others/population_files_HO folder
```


15. Create a merging file that is necessary to merge all populations together
```
# Check if file already exists. If it does, delete it (otherwise we append it too much)
file_to_delete=~/plink/merging_files/populations_file_list.txt
if [ -f "$file_to_delete" ]; then
    rm "$file_to_delete"
    echo "File $file_to_delete deleted."
else
    echo "File $file_to_delete does not exist."
fi

# Set the directory path
directory_path="/home/ubuntu/plink/data_others/extracted_SNPs"

# Create an array to store population file lists
declare -a populations

# Iterate through the files in the directory
for file in "$directory_path"/*.bed; do
    if [ -f "$file" ]; then
        # Extract population name from the filename
        population=$(basename "$file" .bed)

        # Append the file names to the array in the desired order
        populations+=("$population.bed $population.bim $population.fam")
    fi
done

# Print the population file lists, each on a new line
printf '%s\n' "${populations[@]}" > ~/plink/merging_files/populations_file_list.txt

# We add the controls as well!
echo /home/ubuntu/plink/our_samples/extracted_SNPs/unrelated_controls.bed /home/ubuntu/plink/our_samples/extracted_SNPs/unrelated_controls.bim /home/ubuntu/plink/our_samples/extracted_SNPs/unrelated_controls.fam >> ~/plink/merging_files/populations_file_list.txt

# With this file, we can actually merge all of the Human Origins data
```


16. Merge Human Origins files
```
# Go into the directory with all the files
cd ~/plink/data_others/extracted_SNPs

# Perform the merge
plink --bfile ~/plink/our_samples/extracted_SNPs/unrelated_raute --allow-no-sex \
        --merge-list ~/plink/merging_files/populations_file_list.txt \
        --make-bed \
        --out ~/plink/modified_samples/HO_samples
```

# See next script to add the non-HO data from Arciero et al. (2018)






















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
mkdir ~/plink/modified_samples
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
mkdir ~/plink/merging_files; cd ~/plink/mergingfiles #make new directory
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

7. Format the merged data the way you want. I wanted to make sure all samples have the same phenotype and to change the FID from plink to the population names.
```
plink --bfile ~/plink/modified_samples/all_pops --allow-no-sex \
	--pheno ~/plink/subset_text_files/update_all_phenos.txt \
	--update-ids ~/plink/subset_text_files/AllPops_updateFamilyIDs.txt \
	--make-bed \
	--out pheno_ready
```

3. First 



























