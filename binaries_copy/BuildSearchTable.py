# Code taken from https://stackoverflow.com/questions/73984577/read-csv-files-to-pandas-dataframe-from-ftp-with-python-ftplib
import pandas as pd
import ftplib as ftplib
import Bio
import os
import sys
from io import BytesIO as BytesIO
from os.path import join
from pandas.api.types import CategoricalDtype
from ete3 import NCBITaxa, Tree, TreeStyle, BarChartFace 
from Bio import Entrez

%matplotlib inline
ncbi = NCBITaxa()
# If this doesn't work run this database update!!
# ncbi.update_taxonomy_database()

# Set datatypes
dtypes = {'#Complex ac': 'category',
  'Recommended name': 'category',
  'Aliases for complex': 'category',
  'Taxonomy identifier': 'category',
  'Evidence Code': 'string',
  'Experimental evidence': 'string',
  'Go Annotations': 'string',
  'Cross references': 'string',
  'Description': 'string',
  'Complex properties': 'string',
  'Complex assembly': 'string',
  'Ligand': 'string',
  'Disease': 'string',
  'Agonist': 'string',
  'Antagonist': 'string',
  'Comment': 'string',
  'Source': 'string',
  'Expanded participant list': 'string',
  'Identifiers (and stoichiometry) of molecules in complex': 'string'}

# set exclusion list
exclusion_list = {'562.tsv','83333.tsv','2697049.tsv','694009.tsv','1263720.tsv','1235996.tsv','208964.tsv','243277.tsv'}

site = 'ftp.ebi.ac.uk'
### Function to download TSVs
def fetch_files(site, username, password, directory: str = '/pub/databases/intact/complex/current/complextab/', filematch: str = '*.tsv'):
    with ftplib.FTP(site) as ftp:
        # pass the url without protocol
        ftp = ftplib.FTP(site)
        # pass credentials if anonymous access is not allowed
        ftp.login(username, password)
        ftp.cwd(directory)
        list_ = []
        for file_ in ftp.nlst(filematch):
            if file_ in exclusion_list:
                continue
            else:
    #           print(file_) # This works
                flo = BytesIO()
                ftp.retrbinary('RETR ' + file_, flo.write)
                flo.seek(0)
                df = pd.read_csv(flo, na_values=['_'], sep='\t', 
                dtype = dtypes
                # Change column data types (Python Course Lesson 1 - cell 45) 
    #            dtype = {name:'string' for name in Complexes2.select_dtypes('object').columns}
            )
            list_.append(df)
    return pd.concat(list_)
Complexes = fetch_files('ftp.ebi.ac.uk','','')

# Remove stoichiometric data and split accessions in ID column
Complexes['Identifiers (and stoichiometry) of molecules in complex'] = Complexes['Identifiers (and stoichiometry) of molecules in complex'].str.replace('\(.\)','').str.split('|')
# Flatten Taxon ID list and move it to end of duplicate dataframe
flattened_col = pd.DataFrame([(index, value) for (index, values) in Complexes['Identifiers (and stoichiometry) of molecules in complex'].iteritems() for value in values],
                             columns=['index', 'Identifiers (and stoichiometry) of molecules in complex']).set_index('index')
Complexes2 = Complexes.drop('Identifiers (and stoichiometry) of molecules in complex', axis=1).join(flattened_col)
Complexes2.dtypes
Complexes2

# Uncomment the next line to write the whole dataframe to a file
# Complexes2.to_csv('Complexes2.csv')

# Rewrite column dtypes
# dtype = {name:'string' for name in Complexes2.select_dtypes('object').columns},{name:'string' for name in Complexes2.select_dtypes('int64').columns}

# ---> IF POSSIBLE, THIS IS WHERE IT WOULD BE GOOD TO EXECUTE THE ACTUAL SEARCH TO RETRIEVE THE RESULTS TABLE FOR THE NEXT PART (blastp)

# Load in results table - taxid needs to be string to split it later
dtypesa = {'query': 'category',
           'taxid': 'string',
           'orthologs': 'string',
           'identity': 'float64',
           'evalue': 'float64',
           'bitscore': 'float64'}

blastp = pd.read_csv(
    'Diamond_blastp_nr_results.tsv',
    sep="\t",
    dtype=dtypesa
    )

# Drop columns for ease of use
Complexes2 = Complexes2.drop(columns=['Aliases for complex', 'Evidence Code', 'Experimental evidence', 'Go Annotations', 'Cross references', 'Description', 'Complex properties', 'Complex assembly', 'Ligand', 'Disease', 'Agonist', 'Antagonist', 'Comment', 'Source', 'Expanded participant list'])

# Change datatypes to category
Complexes2[Complexes2.columns[:4]] = Complexes2[Complexes2.columns[:4]].astype('category')
blastp[blastp.columns[:3]] = blastp[blastp.columns[:3]].astype('category')
blastp[blastp.columns[1]] = blastp[blastp.columns[1]].astype('string').str.split(';')

#Mergeattempt = Complexes2.merge(
#    blastp,
#    left_on=(['Identifiers (and stoichiometry) of molecules in complex']),
#    right_on=(['query']),
#    how='left'
#)

#20230608 - This is where particular complexes can be subsetted based on name
pattern = 'mcm2|mcm3|mcm4|mcm5|mcm6|mcm7|Cdc45|RPA1|primase|DNA polymerase|RFC1|RFC2|RFC3|RFC4|RFC5|PCNA|Fen1'
Mergeattempt = Complexes2[Complexes2['Recommended name'].str.contains(pattern, case=False, na=False)].merge(
    blastp,
    left_on=(['Identifiers (and stoichiometry) of molecules in complex']),
    right_on=(['query']),
    how='left'
)

# Retrieve taxids
Treetaxa = Mergeattempt['taxid'].explode().unique()

# Write tree
# Need to drop NAs first
Treetaxa = Treetaxa[~pd.isnull(Treetaxa)]

# If this doesn't work because of Keyerrors, find out which ones they are:
# Treetaxa.tolist().index('2461416')

# Remove them!
l = Treetaxa.tolist()
#l.remove('2461416')
#l

# Try again; THIS WORKS!
tree = ncbi.get_topology((l[:10000]),intermediate_nodes=False)
tree.write(format=9,features=["sci_name", "rank"],format_root_node=False,outfile='Treetaxa.nwk')

# Rename tree nodes - May need to comment out or do this after connecting the heatmap
#for node in tree.traverse():
#        node.name = node.sci_name

# Can render the tree inline
# tree.render("%%inline")

from ete3 import ClusterTree
matrix = Mergeattempt[['taxid','orthologs']]


# Make a search command to generate presence absence heatmap (from Martin): Note for this version I removed the > 0 that would make it boolean. This is pretty much a measure of number of hits, which range from 654 to 1
matrix = (Mergeattempt[['#Complex ac', 'taxid']][:500000000].explode('taxid')[['#Complex ac', 'taxid']].value_counts().unstack(0).fillna(0)).reset_index()


matrix.columns.name = None
matrix.set_index('taxid')

#matrix_string = matrix.rename(columns={'taxid': '\#Names'}).head(200000).to_csv(None, sep='\t', index=False)[1:]

# For ClusterTree, tree has to be a newick string (can't be phylo object), matrix has to be a string
matrix_string = matrix.rename(columns={'taxid': '\#Names'}).to_csv(None, sep='\t', index=False)[1:] #matrix

# WORKABLE SUBSET
# Troubleshooting
t = ClusterTree('/Users/of2/Documents/Complex_Portal/ComPred/Eukaryotic_subset.nwk', text_array=matrix_string)

# Lepidoptera
#tree = ncbi.get_topology([7088])
#u = tree.write()
t = ClusterTree(u,text_array=matrix_string)

# Rename tree nodes - May need to comment out or do this after connecting the heatmap
#ncbi.annotate_tree(t, taxid_attr="name")

# Try and get labels - from https://stackoverflow.com/questions/51366712/adding-labels-to-heatmaps-in-ete3
labels = matrix_string.split('\n')[0].split('\t')[1:]
axisface = BarChartFace([0]*len(labels), width=320, height=0, labels=labels, label_fsize=1, max_value=0.5, scale_fsize=6)
axisface.margin_bottom = 30

ts = TreeStyle()
ts.show_leaf_name = True 

ts.aligned_foot.add_face(axisface, column=0)
ts.show_border = False
ts.margin_top = 50
for node in t.traverse():
    node_name = str(join(node.name))
    if node_name:  # Check if node_name is not empty
        node.name = ncbi.get_taxid_translator([node_name])
#ncbi.annotate_tree(t, taxid_attr="name")
t.render("%%inline", 'heatmap', w=1583, units="mm", tree_style=ts)
#################################################
# Full set

t = ClusterTree('Treetaxa.nwk', text_array=matrix_string)
                #matrix.head().to_csv(None, sep='\t')) 
# Try and get labels - from https://stackoverflow.com/questions/51366712/adding-labels-to-heatmaps-in-ete3
labels = matrix_string.split('\n')[0].split('\t')[1:]
axisface = BarChartFace([0]*len(labels), width=535, height=-30, labels=labels, max_value=1, scale_fsize=4)
axisface.margin_bottom = 10
ts = TreeStyle()
ts.show_leaf_name = True
ts.aligned_header.add_face(axisface, column=0)
ts.show_border = False
ts.margin_top = 50
for node in t.traverse():
    node_name = str(join(node.name))
    if node_name:  # Check if node_name is not empty
        node.name = ncbi.get_taxid_translator([node_name])
t.render("%%inline", 'heatmap', tree_style=ts)
t.show("heatmap")
t.show("cluster_cbars")
t.show("cluster_bars")
t.show("cluster_lines")
# Write Merge attempt
# Mergeattempt.to_csv('Mergeattempt.tsv', sep='\t')

# 20230608 - my attempt at filtering complexes to get subset
pattern = 'mcm2|mcm3|mcm4|mcm5|mcm6|mcm7|Cdc45|RPA1|primase|DNA polymerase|RFC1|RFC2|RFC3|RFC4|RFC5|PCNA|Fen1'
Complexes2[Complexes2['Recommended name'].str.contains(pattern, case=False, na=False)]

theirtext ="""
#Names\tcol1\tcol2\tcol3\tcol4\tcol5\tcol6\tcol7
A\t-1.23\t-0.81\t1.79\t0.78\t-0.42\t-0.69\t0.58
B\t-1.76\t-0.94\t1.16\t0.36\t0.41\t-0.35\t1.12
C\t-2.19\t0.13\t0.65\t-0.51\t0.52\t1.04\t0.36
D\t-1.22\t-0.98\t0.79\t-0.76\t-0.29\t1.54\t0.93
E\t-1.47\t-0.83\t0.85\t0.07\t-0.81\t1.53\t0.65
F\t-1.04\t-1.11\t0.87\t-0.14\t-0.80\t1.74\t0.48
G\t-1.57\t-1.17\t1.29\t0.23\t-0.20\t1.17\t0.26
H\t-1.53\t-1.25\t0.59\t-0.30\t0.32\t1.41\t0.77
"""

mytext="""
#Names\tCPX-2635\tCPX-15\tCPX-259\tCPX-365\tCPX-33\tCPX-1291\tCPX-540\tCPX-6\tCPX-543\tCPX-2187\tCPX-1154
100035\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t1.0
1000413\t17.0\t17.0\t17.0\t17.0\t17.0\t17.0\t17.0\t17.0\t17.0\t17.0\t0.0
1001064\t2.0\t2.0\t2.0\t2.0\t2.0\t2.0\t2.0\t2.0\t2.0\t2.0\t1.0
1001832\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t6.0\t1.0
1001833\t7.0\t7.0\t7.0\t7.0\t7.0\t7.0\t7.0\t7.0\t7.0\t7.0\t1.0
"""

# Extract the 'Go Annotations' column
go_annotations = Complexes['Go Annotations']
# Split the multiple entries delimited by a pipe symbol
split_annotations = go_annotations.str.split('|')
# Flatten the list of annotations
flat_annotations = [annotation for sublist in split_annotations.dropna() for annotation in sublist]
# Count the frequency of each annotation
annotation_counts = pd.Series(flat_annotations).value_counts()
# Print the abundance of each gene ontology annotation
print(annotation_counts)

ligase_entries = annotation_counts[annotation_counts.index.str.contains('ligase', case=False)]
print(ligase_entries)



for node in t.traverse():
    node_name = str(join(node.name))
    if node_name:  # Check if node_name is not empty
        print(ncbi.get_taxid_translator([node_name]))o

Tetrahymena thermophila instead of Oxytricha - 5911 instead of 1172189
Theileria parva muguga instead of Theileria parva - 333668 instead of 5875
Allomyces macrogynus ATCC 38327 instead of Allomyces macrogynus - 578462 instead of 28583

Spizellomyces punctatus DAOM BR117 instead of Spizellomyces punctatus - 645134 instead of 109760
Ustilago maydis 521 (smut fungi) instead of Ustilago maydis - 237631 instead of 5270
Heterostelium album PN500 instead of Heterostelium pallidum - 670386 instead of 13642
Capsaspora owczarzaki ATCC 30864 instead of Capsaspora owczarzaki - 595528 instead of 192875

Monosiga brevicollis MX1 instead of Monosiga brevicollis - 431895 instead of 81824

Guillardia theta CCMP2712 instead of Guillardia theta - 905079 instead of 55529
Emiliania huxleyi CCMP1516 instead of Emiliania huxleyi - 280463 instead of 2903
instead of Encephalitozoon cuniculi

Bigelowiella natans CCMP2755 instead of Bigelowiella natans - 753081 instead of 227086