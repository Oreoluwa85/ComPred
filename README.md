# ComPred - Prediction of eukaryotic protein complexomes by homology searching

### Description
ComPred is a series of scripts to produce a taxonomy tree adjacent heatmap, corresponding to presence or absence of any protein component of protein complexes identified in the [EBI Complex Portal](https://www.ebi.ac.uk/complexportal/home) (recent citation can be found at [Meldal et al, 2021](https://www.doi.org/10.1093/nar/gky1001)). 

[//]: # (What is this idea based upon - Koonin bifurcation theory, inferring function from complex conservation)
[//]: # (, the domain of life comprising fungi, plants, ourselves and the rest of the animals as well as the lesser-mentioned protists)

Despite known occurences of non-vertical evolution ([Gabaldon 2020](https://doi.org/10.1111/nyas.14471)), comparative ancestral genome analysis remains a valuable tool in the ongoing characterisation of eukaryotes. The ComPred toolset uses **DIAMOND BLASTp** to search any sequence database (here, in developmental stages, the non-redundant NCBI database) for homologs to FASTA sequences corresponding to each component of each complex in the EBI Complex Portal.

### Contents

**BuildSearchTable.py**: The python only version downloads the latest release of complex portal data as TSVs files, then concatenates them into a single table. A FASTA file corresponding to protein sequences for each protein in this table is used to generate a search a protein sequence database using DIAMOND BLASTp. The file then takes a newick file of NCBI taxids and aligns with presence or absence values of homologs to each Comple Portal complex, where a presence of any component corresponds to a complex presence and an absence of all components corresponds to a complex absence. The image can be viewed inline, or output as a png using ete3.tree.show functionality.

**diamond_blastp_query.sh**: The protein sequence search query code. This default searches a local copy of the non-redundant NCBI database, returning a single hit for each query (each protein from the Complex Portal) for each taxon. The resulting table is a large, parsable TSV to be queried with **BuildSearchTable.py**.

### Old Scripts
* Earlier scripts required python3, SQLite3 and Perl and have been rewritten and superceded by a python only version. The old scripts should still work, but are not recommended. Instructions are at the bottom of the page.

The old scripts are to be called as follows:

**Complex_Portal_Retrieval**.sh: Download latest release of complex portal data as TSVs and corresponding FASTA files, then concatenate each into a single query file. The FASTA file is for querying the database and the TSV (which has complex and protein accession numbers) is for indexing the results of the sequence search.

**diamond_blastp_query.sh**: This is the diamond protein search script that searches a DIAMOND database for the FASTA sequences in the concatenated file from the previous step.

**Concat_results**: This script concatenates and removes duplicate hits of each query FASTA for each taxon from the resulting DIAMOND BLAST or HMMER3 results table.

**make_compdb.sh**: This script creates a SQLite3 database, with the concatenated Complex_Portal_Retrieval.sh indexing table and the de-duplicated DBLAST/H3 results table.

**Get_taxid_newick.py**: This script creates a newick tree from a text file contaning a list of NCBI taxon IDs.

**make_resultsgrid.r**: This script creates the resultsgrid.csv, a matrix of presence or absence of protein complexes.

[//]: # (NEED TO:)
[//]: # (	MAKE AN R SCRIPT THAT TAKES IN THE OUTPUT SQLITE3 DB FILENAME FROM MAKE_COMPDB.SH)
[//]: # (	MAKE A SCRIPT THAT FINDS THE CORRESPONDING ALIGNMENT FROM EACH UNIPROT ACCESSION AND TO MAKE AN ALIGNMENT, WHICH CAN BE USED TO MAKE A PROFILE)
[//]: # (	MAKE A SCRIPT THAT QUERIES EACH PROTEIN IN EACH PROTEOME ATTACHED TO EACH GENOME AGAINST EACH PROTEIN'S HMM)
