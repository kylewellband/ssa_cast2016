#!/bin/bash

mkdir ./02_Results/00_Basic_Stats
pops=`cut -f1 -d " " ./01_Data/CAST_2016_50K.ped | uniq`

for p in $pops;
do
	grep -e $p ./01_Data/CAST_2016_50K.ped > ./01_Data/${p}.ped;
	/Applications/plink_mac/plink --ped ./01_Data/${p}.ped --map ./01_Data/CAST_2016_50K.map --chr-set 29 --not-chr 8,29 --hardy midp --out ./02_Results/00_Basic_Stats/${p}.tmp
	grep -e "ALL" ./02_Results/00_Basic_Stats/${p}.tmp.hwe > ./02_Results/00_Basic_Stats/${p}.hwe
	rm ./02_Results/00_Basic_Stats/${p}.tmp.hwe
	/Applications/plink_mac/plink --ped ./01_Data/${p}.ped --map ./01_Data/CAST_2016_50K.map --chr-set 29 --not-chr 8,29 --het --ibc --out ./02_Results/00_Basic_Stats/${p}
	rm ./01_Data/${p}.ped
done
