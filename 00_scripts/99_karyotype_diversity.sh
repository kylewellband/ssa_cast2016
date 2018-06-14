#!/bin/bash

## 
## by: Kyle Wellband
## Script to generate diversity within and divergence among karyotypes using PLINK


cd ~/Desktop/CAST_2016_50K_analysis

scripts_dir="./00_scripts/"
data_dir="./01_data/"
results_dir="./02_results/"
misc_dir="./05_misc/"

# Split dataset by inversion/fusion gt
for gt in AA AB BB; do
    grep -e $gt ${data_dir}/CAST_2016_50K_inversion_status.txt | cut -f1,2 -d " " > ${data_dir}/indsub.tmp
    /Applications/plink_mac/plink --bfile ${data_dir}/chr8_29 --chr-set 29 --keep ${data_dir}/indsub.tmp --make-bed --hardy --freq --r2 square yes-really --out ${data_dir}/chr8_29_$gt
    rm ${data_dir}/indsub.tmp
done

# Overall FST
/Applications/plink_mac/plink --bfile ${data_dir}/chr8_29 --chr-set 29 --fst --within ${data_dir}/CAST_2016_50K_inversion_status.txt --out ${data_dir}/CAST_2016_50K_inversion

# pairwise FST
/Applications/plink_mac/plink --bfile ${data_dir}/chr8_29 --chr-set 29 --fst --within ${data_dir}/CAST_2016_50K_inversion_status.txt --remove-cluster-names AB --out ${data_dir}/CAST_2016_50K_inversion_AA_BB
/Applications/plink_mac/plink --bfile ${data_dir}/chr8_29 --chr-set 29 --fst --within ${data_dir}/CAST_2016_50K_inversion_status.txt --remove-cluster-names AA --out ${data_dir}/CAST_2016_50K_inversion_AB_BB
/Applications/plink_mac/plink --bfile ${data_dir}/chr8_29 --chr-set 29 --fst --within ${data_dir}/CAST_2016_50K_inversion_status.txt --remove-cluster-names BB --out ${data_dir}/CAST_2016_50K_inversion_AA_AB

