# Merging SNP array data with PLINK
Here, I merge data from different types of arrays, e.g. different Illumina arrays and the Human Origins array, to produce a large dataset. Almost all data I use is public and can be obtained via https://reich.hms.harvard.edu/datasets. Some datasets I use for my analyses are not public and therefore I do not show how to download any datasets, but only call them by first author and year. I also run IBIS on all of these datasets to make sure that the individuals included are minimally related or at least in the same way as the PI_HAT threshold I set for my own samples. 

1. Create directory for data from other labs.
```
mkdir ~/plink/data_others
# Go into directory and download all data that you want to add using for instance 'wget <link>' if public. If private, use way specified by authors/committee.
```

2. All datasets need to be filtered in the same way
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

3. We have to check if any of these samples still contain any related individuals according to the threshold we set. Add map to all the bim files in the directory so that we can run IBIS.
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

4. Run IBIS on each of these datasets.
```
# List of file names (without extensions)
files=("clean_Wang_2021" "clean_Lazaridis_2014" "clean_Jeong_2019" "clean_Skoglund_2016" "clean_Changmai_2022" "clean_Nakatsuka_2017" "clean_Agta_2023" "clean_Arciero_2018")

# Loop through the file names
for file in "${files[@]}"; do
    # Define the input file names
    bed_file="${file}.bed"
    bim_file="${file}.bim"
    fam_file="${file}.fam"

    # Run the IBIS program with the specified parameters
    ~/ibis/ibis "$bed_file" "$bim_file" "$fam_file" -min_l 7 -mt 500 -2 -mt2 500 -er .004 -printCoef -f ~/out

    # Optionally, you can add more commands or processing here for each file
done
```
