#!/bin/bash
read -p "This script downloads ComplexTab versions of complexome proteins from the EBI Complex Portal, resulting in a single FASTA file of all components for DIAMOND, HMMER3 or MMSeq2 search. If these files have been downloaded already, cancel this script now with Ctrl+c"

## Download ComplexTab versions of complexomes
cat TaxonLists/Complex_Portal_Taxonomy_Accessions.md| while read p;do
#          wget "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/$p.tsv"
# Curl if Wget doesn't work
	curl --output $p.tsv "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/$p.tsv"
  done



## Store Uniprot accessions from each Complex Portal tsv file as a txt file
#/\(Q\|P\|B\|O\)[0-9]\S\{3\}[0-9]
#/\([A-Z]\)[0-9]\S\{3,5\}[0-9](.)
for i in *.tsv; do
        grep -o '[A-Z][0-9]\S\{3,5\}[0-9](.)' "$i" >> "${i/%.tsv/.txt}";done

## Remove numerator in the brackets
#for i in *.txt; do; sed -i "" 's/([0-9])$//g' "$i"; done

for i in *.txt; do
        sed -i '' -e 's/([0-9])$//g' "$i"
done

## Retrieve FASTA files from uniprot IDs corresponding to each taxa with a complexome (for BLASTING other taxa)
for i in *.txt; do
  while read p;do
    sort $i | uniq -u
    curl -s "https://rest.uniprot.org/uniprotkb/$p.fasta" >> "$i.fasta"
  done < "$i"
done

## Rename files in Folder
mv 2697049.txt.fasta SARS-CoV-2.txt.fasta
mv 1263720.txt.fasta BetaCoV1.txt.fasta
mv 1235996.txt.fasta HuBetaCoV2c.txt.fasta
mv 694009.txt.fasta SARS-CoV.txt.fasta
mv 559292.txt.fasta S_cerevisiae.txt.fasta
mv 284812.txt.fasta S_pombe.txt.fasta
mv 243277.txt.fasta V_cholerae.txt.fasta
mv 208964.txt.fasta P_aeruginosa.txt.fasta
mv 83333.txt.fasta E_coli.txt.fasta
mv 10116.txt.fasta R_norvegicus.txt.fasta
mv 10090.txt.fasta M_musculus.txt.fasta
mv 9986.txt.fasta O_cuniculus.txt.fasta
mv 9940.txt.fasta O_aries.txt.fasta
mv 9913.txt.fasta B_taurus.txt.fasta
mv 9823.txt.fasta S_scrofa.txt.fasta
mv 9615.txt.fasta C_lupus.txt.fasta
mv 9606.txt.fasta H_sapiens.txt.fasta
mv 9031.txt.fasta G_gallus.txt.fasta
mv 8732.txt.fasta C_durissus.txt.fasta
mv 8355.txt.fasta X_laevis.txt.fasta
mv 7955.txt.fasta D_rerio.txt.fasta

## Concatenate fastas
cat *.txt.fasta >> Allseqs.fasta

## Remove error messages; lines starting with "Error" or "The"
#sed -i 's/\(Error\|The\).*$//' Allseqs.fasta
sed -i '' -e 's/\(.\)>/\1\n>/' Allseqs.fasta
sed -i '' -e 's/\(Error\).*$//' Allseqs.fasta
sed -i '' -e 's/\(The\).*$//' Allseqs.fasta
sed -i '' -e 's/\(Resource\).*$//' Allseqs.fasta

## Join sequences to single line
sed -e 's/\(^>.*$\)/#\1#/' Allseqs.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' >> Allseqs_concat.fasta

## Remove duplicate sequences
sed -e '/^>/s/$/@/' -e 's/^>/#/' Allseqs_concat.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t $'\t' -f -k 2,2  | sed -e 's/^/>/' -e 's/\t/\n/' > Allseqs_concat_uniq.fasta
echo "Query FASTAs ready"

## Concatenate Complex Portal tables
cat *.tsv > Complex_Portal_Table_raw.tsv
sed -i '' -e '2,$s/#.*$//' Complex_Portal_Table_raw.tsv
sed -i '' -e '/^$/d' Complex_Portal_Table_raw.tsv
sed -i '' -e '1s/#//' Complex_Portal_Table_raw.tsv
echo "Query Complex Tables ready. Unnesting accessions in Complex Table"
./Amend_Complex_Table.R
echo "Done"
echo "if parsing Eggnog Results us the following:\n./EggnogResultParse.sh Complex_portal_EggnogResults2.tsv > emapper_parse.tsv"
