from ete3 import NCBITaxa
from Bio import Entrez
import os
import sys
ncbi = NCBITaxa()
# If this doesn't work run this database update!!
# ncbi.update_taxonomy_database()

try:
	if len(sys.argv) > 2:
		filename = sys.argv[1]
		filename2 = sys.argv[2]
		print('Parsing files from ' + filename + ' into ' + filename2)
		if not filename:
			print('Usage python ./Get_taxid_newick.py <taxon.list> <treefile.nwk>')
			exit()
		elif not filename2: 
			print('Usage python ./Get_taxid_newick.py <taxon.list> <treefile.nwk>')
			exit()
		else:
			print('thanks')
	else:
		print('Usage python ./Get_taxid_newick.py <taxon.list> <treefile.nwk>')
		exit()
except IndexError:
		print('Usage python ./Get_taxid_newick.py <taxon.list> <treefile.nwk>')
		exit()
 
try: 
	with open(filename) as f:
		lines = f.read().splitlines()
		tree = ncbi.get_topology((lines),intermediate_nodes=False)
		tree.write(format=9,features=["sci_name", "rank"],format_root_node=False,outfile=filename2)
except FileNotFoundError:
		print('Usage python ./Get_taxid_newick.py <taxon.list> <treefile.nwk> - taxon list not found')
