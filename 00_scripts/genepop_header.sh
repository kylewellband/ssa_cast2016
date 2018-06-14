#!/bin/bash
# script to add/edit header of a genepop file

usage="Usage: $0 <\"header text\" [string]> <filename>"

if [ $# -ne 2 ]; then 
    echo
    echo "Error: incorrect number of arguments"
    echo
    echo $usage
    echo
    exit 1 
fi

cat $2 | tail -n+2 > $2.tmp

echo '"'$1'"' | cat - $2.tmp > $2

rm $2.tmp