#!/bin/bash
filename=$1
filename2=$2
filename3=$3
if [ -z "$filename" ]
then
	echo "usage = ./make_compdb.sh <ComplexTable> <Results/BLAST/HMMER Table> <DB filename>"
	exit
elif [ -z "$filename2" ]
then
	echo "usage = ./make_compdb.sh <ComplexTable> <Results/BLAST/HMMER Table> <DB filename>"
	exit
else
	echo "Building Sqlite database. First table is complex portal components. Second table is search results"
	echo -e ".mode tabs\n.import '$1' Complex_portal_components\n.import '$2' Complex_portal_Results" > temp_make_complex.sql
fi

sqlite3 $filename3 < temp_make_complex.sql
rm temp_make_complex.sql
#sqlite3 ComplexPortalSearch2.db < make_complex_database.sql
#sqlite3 ComplexPortalSearch3.db .mode tabs
#sqlite3 ComplexPortalSearch3.db .import '$1' Complex_portal_components
#sqlite3 ComplexPortalSearch3.db .import '$2' Complex_portal_Results
