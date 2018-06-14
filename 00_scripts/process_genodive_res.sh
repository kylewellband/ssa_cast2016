#!/bin/bash
# script to process genodive results to create files of individuals that have been
# excluded using permutation <$1.excluded> and assignment results <$1.results>

#sed -e '1,/Ind/ d' -e '/^$/,$ d' results_training1_500.gdv > test

#sed -e '/Ind/,$ p' -e '/^$/,$ d' results_training1_500.gdv > test
if [ $# -ne 1 ]; then 
    echo
    echo "Error: No genodive results file provided"
    echo
    exit 1
fi

out=$(echo $1 | sed 's/\.gdv//')

cat $1 | sed -n '/Ind/,$ p' | sed '/^$/,$ d' > $out.excluded

cat $1 | sed '1,/Ind/ d' | sed -n '/Ind/,$ p' | sed '/^$/,$ d' > $out.results

