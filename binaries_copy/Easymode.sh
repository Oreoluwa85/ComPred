#!/bin/bash
Complex_Portal_Retrieval.sh
diamond_blastp_query.sh
Concat_results.sh SupfamResults.tsv
make_compdb.sh Complex_Portal_Table.tsv  SupfamResults.tsv ComplexPortalSearch_20230303
python Get_taxid_newick.py <MarkGSheet.list> <Automated.nwk>
./make_resultsgrid.r
