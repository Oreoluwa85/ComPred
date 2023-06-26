#!/bin/bash
read -p "This script performs a Diamond Search with a FASTA file of EBI Complex Portal components ran via Complex_Portal_Retrieval.sh. If you are not on the farm, if you don't have Diamond installed or if you don't want to run this, cancel this script now with Ctrl+c. If you already have the results from running this, cancel and run Concatenate_results.sh "

echo "this job is running on:"
hostname
/software/team301/./diamond \
      blastp \
      --threads 8 \
      --db /lustre/scratch123/tol/teams/blaxter/users/rc28/2021/databases/uniprot_2021_06/reference_proteomes.dmnd \
      --query /lustre/scratch123/tol/teams/blaxter/users/of2/Allseqs_acc.fa \
      --outfmt 6 qseqid staxids sseqid pident evalue bitscore \
      --evalue 1.0  --max-target-seqs 3000 --max-hsps 1 \
      --block-size 6 \
      --out /lustre/scratch123/tol/teams/blaxter/users/of2/Diamond_blastp_results.txt
