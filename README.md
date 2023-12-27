# Merging non-Human Origins data

Here, we continue with merging our data. We want to add the data from Arciero et al (2018) with our other data. However, this was genotyped on an Illumina array, hence merging will look slightly differently.

1. Sort SNPs 
```
echo "Start step 1"

# Move into directory
cd ~/plink/data_others/population_files_ARC/

# Make new directory
mkdir ~/plink/data_others/population_files_ARC/sorted_SNPs

# List populations
find ~/plink/data_others/population_files_ARC/ -type f -name '*.bim' | awk -F/ '{print $NF}' | sort -u > files

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

done < files

# Move all files into new directory
mv ./*_sorted_SNPs.txt ~/plink/data_others/population_files_ARC/sorted_SNPs
echo "sorted_SNPs.txt moved into sorted_SNPs directory"

echo "End step 1"
```


2. Combine sorted SNPs
```
echo "Start step 2"

cd ~/plink/data_others/population_files_ARC/sorted_SNPs

# Combine all files into one file
cat ./*_sorted_SNPs.txt > combined_snps.txt
echo "combined_snps.txt created"

# Sort this file
sort combined_snps.txt -o combined_snps.txt
echo "combined_snps.txt ordered"

# Find all unique files and sort them by frequency
uniq -c combined_snps.txt > frequency_snps.txt
echo "frequency_snps.txt created"

# Find the number of populations so that we can select only the SNPs which occur in all pops
awk '{if ($1 > max) max = $1} END {print max}' frequency_snps.txt

# This should be 12 (as of Tue 19 December 2023).

# Select only the SNPs which occur in all pops
awk '$1 == 24 {print $2}' frequency_snps.txt > overlapping_snps.txt
echo "overlapping_snps.txt created"

# This gives: 158011 (Tue 19 Dec 2023) SNPs which they occur in EACH single population!
wc -l ~/plink/data_others/population_files_ARC/sorted_SNPs/overlapping_snps.txt

echo "End step 2"
```


3. Extract those SNPs
```
echo "Start step 3"

# Working directory
cd ~/plink/data_others/population_files_ARC/

# Make new directory
mkdir ~/plink/data_others/population_files_ARC/extracted_SNPs

# Use 'files' from step 1.

echo "Start loop to extract common SNPs"

# Perform IBIS for each population
while IFS= read -r file; do

file_prefix=$(basename "$file" | cut -d. -f1)
plink2 --bfile ./$file_prefix \
		--extract ~/plink/data_others/population_files_ARC/sorted_SNPs/overlapping_snps.txt \
		--make-bed \
		--out ~/plink/data_others/population_files_ARC/extracted_SNPs/$file_prefix

done < files

echo "End loop"

# Remove files in the ~/plink/data_others/population_files_ARC folder
```


4. Create merge file to merge Arciero populations again
```
# Check if file already exists. If it does, delete it (otherwise we append it too much)
file_to_delete=~/plink/merging_files/populations_arc_file_list.txt
if [ -f "$file_to_delete" ]; then
    rm "$file_to_delete"
    echo "File $file_to_delete deleted."
else
    echo "File $file_to_delete does not exist."
fi

# Set the directory path
directory_path="/home/ubuntu/plink/data_others/population_files_ARC/extracted_SNPs"

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
printf '%s\n' "${populations[@]}" > ~/plink/merging_files/populations_arc_file_list.txt

head -n -1 ~/plink/merging_files/populations_arc_file_list.txt > ~/plink/merging_files/populations_arc_sub_list.txt
rm ~/plink/merging_files/populations_arc_file_list.txt
```


5. Merge Arciero files
```
# Go into the directory with all the files
cd ~/plink/data_others/population_files_ARC/extracted_SNPs/

# Perform the merge
plink --bfile Tingri --allow-no-sex \
        --merge-list ~/plink/merging_files/populations_arc_sub_list.txt \
        --make-bed \
        --out ~/plink/modified_samples/ARC_samples
```


6. Find overlapping SNPs between HO data and Arciero data
```
cd ~/plink/merging_files

# Find and extract the IDs common in both files
cat ~/plink/data_others/sorted_SNPs/overlapping_snps.txt ~/plink/data_others/population_files_ARC/sorted_SNPs/overlapping_snps.txt > combined_arc.txt
sort combined_arc.txt -o combined_arc.txt
uniq -c combined_arc.txt > frequency_arc.txt
awk '$1 == 2 {print $2}' frequency_arc.txt > overlapping_snps_arc.txt
```


7. Merge HO and Arciero data
```
# Go into the directory with all the files
cd ~/plink/data_others/population_files_ARC/extracted_SNPs

# Perform the merge
plink --bfile ~/plink/modified_samples/HO_samples --allow-no-sex \
        --bmerge ~/plink/modified_samples/ARC_samples.bed ~/plink/modified_samples/ARC_samples.bim ~/plink/modified_samples/ARC_samples.fam \
	--extract ~/plink/merging_files/overlapping_snps_arc.txt \
        --make-bed \
        --out ~/plink/merging_files/nonflipped

# The merge doesn't work:
# "Error: 28978 variants with 3+ alleles present.
# * If you believe this is due to strand inconsistency, try --flip with /home/ubuntu/plink/merging_files/nonflipped-merge.missnp. (Warning: if this seems to work,
# strand errors involving SNPs with A/T or C/G alleles probably remain in your data.  If LD between nearby SNPs is high, --flip-scan should detect them.)
# * If you are dealing with genuine multiallelic variants, we recommend exporting that subset of the data to VCF (via e.g. '--recode vcf'), merging with
# another tool/script, and then importing the result; PLINK is not yet suited to handling them."

# I believe it's due to strand inconsistency, hence I flip the samples.

# Flip samples where necessary
plink --bfile ~/plink/modified_samples/ARC_samples \
	--flip ~/plink/merging_files/nonflipped-merge.missnp \
	--make-bed \
	--out ~/plink/merging_files/flipped

# Perform merge with flipped samples
plink --bfile ~/plink/modified_samples/HO_samples --allow-no-sex \
        --bmerge ~/plink/merging_files/flipped.bed ~/plink/merging_files/flipped.bim ~/plink/merging_files/flipped.fam \
	--extract ~/plink/merging_files/overlapping_snps_arc.txt \
	--make-bed \
	--out ~/plink/modified_samples/all_pops
```































  




