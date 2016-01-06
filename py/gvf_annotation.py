import pandas as pd
import sys
from bx.intervals.intersection import Intersecter, Interval
from IPython import embed
import re

# usage: python gvf_annotation.py ch_grch38_uniq_gene_id.gvf duplication_dbVar.gvf 5 RefSeq_Gene_2  
# usage: python gvf_annotation.py 
# command line fields: 1) gvf file with annotations 2) file to be annotated 3) max number of gene matches 4) output file 
print(sys.argv[1])
print(sys.argv[2])

# check if max_matches is defined otherwise set it to 10
try:
     print(sys.argv[3])
     max_matches = int(sys.argv[3])
except IndexError:
     max_matches =10

# check if the output file is set
try:
     print(sys.argv[4])
     output_file = sys.argv[4]
except IndexError:
     output_file = "RefSeq_Gene.gvf"

f1 = pd.read_csv(sys.argv[1], skiprows=1, header=None, sep="\t")

# start and end are 1-based
# fields three and six are set to default value '-'
# field 8 is set to default value '.' to meet compatibility with GFF3
f1.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','attributes']


f2 = pd.read_csv(sys.argv[2], header=None,sep="\t")
f2.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand','phase', 'attributes']

tree_dict = {}

tree = Intersecter()

# first input file - gvf file containing the annotation
for chrom in range(1, 23):
    tree_dict['chr' + str(chrom)] = Intersecter()

tree_dict['chrX'] = Intersecter()
tree_dict['chrY'] = Intersecter()

gene_id_re = re.compile("(GeneID=[0-9]+)")
gene_symbol = re.compile("(GeneSymbol=[A-z0-9]+)")


for i, j in f1.iterrows():
    geneid = gene_id_re.search(j[-1]).group(0)
    geneid_new = geneid.replace("GeneID=","")
    genesymbol = gene_symbol.search(j[-1]).group(0)
    genesymbol_new = genesymbol.replace("GeneSymbol=","")
    tree_dict[j['chr']].add_interval(Interval(int(j.start), int(j.end), 
            value={'geneid_new': geneid_new, 'genesymbol_new': genesymbol_new}))

# open file
f = open(output_file,'w')

# second input file - gvf file to be annotated
for i2, j2 in f2.iterrows():
    hit = tree_dict[str(j2['chr'])].find(j2.start, j2.end)
    if len(hit) > 0:
        new_attr = j2[-1]
        new_attr += ';Match='
        
        # limit the number of assigned genes
        if len(hit) > max_matches: 
            hit = hit[0:max_matches]
        else: pass
        
        cnt = 0
        for hiti in hit:
            # check if it is a exact or partial match
            if (hiti.start == j2.start) & (hiti.end == j2.end):
                    exact_match = 'exact'
            else:
                    exact_match = 'partial'
            
            # delimiter for the output
            cnt += 1
            if (cnt < len(hit)):
                    delim = "|"
            else:
                    delim = ""

            # add information to the attribute field
            new_attr += (str(hiti.start) + ':' + str(hiti.end)+ ':' +
                    exact_match + ":" + hiti.value['geneid_new'] + ':' +
                    hiti.value['genesymbol_new'] + delim)

        # write into file
        f.write("\t".join([str(z) for z in j2[0:-1]]) + "\t" + new_attr + "\n")

# close file        
f.close()
