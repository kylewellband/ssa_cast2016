#!/bin/bash

## Code to document CAST 2016 sub basin analyses
## by: Kyle Wellband

# *** Setup workspace ***

cd ~/Desktop/CAST_2016_50K_analysis

scripts_dir="./00_scripts/"
data_dir="./01_data/"
results_dir="./02_results/"
misc_dir="./05_misc/"


# Convert .ped and .map to binary format for faster processing and initial filtering to remove fixed variants (MAF < 0.001) and variants missing in >10% of individuals and individuals with >10% missing data
/Applications/plink_mac/plink --file ${data_dir}/CAST_2016_50K_all --make-bed --recode --chr-set 29 --maf 0.01 --geno 0.1 --mind 0.1 --out ${data_dir}/CAST_2016_50K


# *** Run outlier tests ***

# Bayescan
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K --chr-set 29 --freq --family --out ${data_dir}/CAST_2016_50K
${scripts_dir}/popfreq2bayescan.R ${data_dir}/CAST_2016_50K.frq.strat
mkdir ${results_dir}/07_outliers/bayescan
/Applications/BayeScan2.1/binaries/BayeScan2.1_macos64bits ${data_dir}/CAST_2016_50K.bayescan -od ${results_dir}/07_outliers/bayescan -o pr_odds_10 -j 2 -pr_odds 10
/Applications/BayeScan2.1/binaries/BayeScan2.1_macos64bits ${data_dir}/CAST_2016_50K.bayescan -od ${results_dir}/07_outliers/bayescan -o pr_odds_10000 -j 2 -pr_odds 10000
${scripts_dir}/bayescan_plotting.R

# baypass
java -Xmx6G -jar /Applications/PGDSpider_2.1.1.2/PGDSpider2-cli.jar -inputfile ${data_dir}/CAST_2016_50K.ped -outputfile ${data_dir}/CAST_2016_50K.bayenv -spid ${misc_dir}/ped2bayenv.spid
${scripts_dir}/bayenv2baypass.R
mkdir ${results_dir}/07_outliers/baypass/
/Applications/baypass_2.1/sources/g_baypass -npop 15 -gfile ${data_dir}/CAST_2016_50K.baypass -outprefix ${results_dir}/07_outliers/baypass/baypass_out

# Perform LD pruning for pop gen analysis
# Based on VIF > 2
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K --chr-set 29 --exclude ${data_dir}/CAST_2016_50K_outliers_to_remove_for_neutral.txt --indep 50 5 2 --out ${data_dir}/CAST_2016_50K_neutral

/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K --chr-set 29 --extract ${data_dir}/CAST_2016_50K_neutral.prune.in --make-bed --recode --out ${data_dir}/CAST_2016_50K_neutral
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K --chr-set 29 --extract ${data_dir}/CAST_2016_50K_bayenv.prune.in --recode --out ${data_dir}/CAST_2016_50K_bayenv
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K_neutral --chr-set 29 --recode A --out ${data_dir}/CAST_2016_50K_neutral

# LDNe
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K_neutral --chr-set 29 --not-chr 8,29 --recode --out ${data_dir}/CAST_2016_50K_neutral_Ne
java -jar /Applications/PGDSpider_2.1.1.2/PGDSpider2-cli.jar -inputfile ${data_dir}/CAST_2016_50K_neutral_Ne.ped -outputfile ${data_dir}/CAST_2016_50K_neutral_Ne.txt -spid ${misc_dir}/ped2genepop.spid
${scripts_dir}/genepop_header.sh "CAST_2016_50K_neutral_Ne_genepop" ${data_dir}/CAST_2016_50K_neutral_Ne.txt

# HWE
${scripts_dir}/01_hwe.sh
${scripts_dir}/01_hwe_inbreeding_rpt.R

# FST
java -jar /Applications/PGDSpider_2.1.1.2/PGDSpider2-cli.jar -inputfile ${data_dir}/CAST_2016_50K_neutral.ped -outputfile ${data_dir}/CAST_2016_50K_neutral.arp -spid ${misc_dir}/ped2arl.spid
/Applications/arlecore_macosx/arlecoremac_64bit ${data_dir}/CAST_2016_50K_neutral.arp ${misc_dir}/arl_run_FST.ars
mv ${data_dir}/CAST_2016_50K_neutral.res/ ${results_dir}/CAST_2016_50K_neutral.res/
${scripts_dir}02_plot_pw_fst.R ${results_dir}/CAST_2016_50K_neutral.res/CAST_2016_50K_neutral.xml ${results_dir}02_

#AMOVA
cat ${data_dir}/CAST_2016_50K_neutral.arp ${misc_dir}/AMOVA_structure.txt > ${data_dir}/CAST_2016_50K_neutral_amova.arp
/Applications/arlecore_macosx/arlecoremac_64bit ${data_dir}/CAST_2016_50K_neutral_amova.arp ${misc_dir}/arl_run_AMOVA.ars
mv ${data_dir}/CAST_2016_50K_neutral_amova.res/ ${results_dir}/CAST_2016_50K_neutral_amova.res/

# admixture
mkdir ${results_dir}/05_admixture
cd ${results_dir}/05_admixture
for K in $(seq 1 10); do /Applications/admixture_macosx-1.3.0/admixture ../../${data_dir}/CAST_2016_50K_neutral.ped --cv ${K} -j1 | tee log${K}.out;done
grep -h CV log*.out > CV_error.txt
cd ../../

# PCA on LD filtered data
/Applications/plink_mac/plink --bfile ${data_dir}/CAST_2016_50K_neutral --chr-set 29 --pca --out ${results_dir}/CAST_2016_50K_neutral
${scripts_dir}/pca_plotting.R ${results_dir}/CAST_2016_50K_neutral

