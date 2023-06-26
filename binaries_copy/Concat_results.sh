#!/bin/bash
read -p "This script retrieves unique hits for each query for each taxon. If you've already ran this, cancel with Ctrl+c. Takes a tsv"
filename=$1
if [ -z "$filename" ]
then
	echo "usage = ./Concat_results.sh <DIAMONDResults.tsv>"
	exit
else
	echo "Removing duplicate taxa specifications"
fi

awk -F"\t" '!uniq[$1, $2]++' $filename
sed -i '1s/^/query\ttaxid\torthologs\tidentity\tevalue\tbitscore\n/' $filename
