import numpy as np
import math
import HTSeq

filepath = "/home/ubuntu/Yan/estd211_Campbell_et_al_2014.variant_call.remap.germline.tab"

def VariantCallTabReader(filepath):
    """
    1. This function is to read from the merged file to obatin the genomic
    regions of each variant call. 
    2. Variant type is assumed to be the same (i.e only deletion at this
    stage). 
    3. No fussiness infomation are considered, i.e only inner
    start and inner end are used.

    """
    infile = open(filepath, 'r')
    # variant_interval is a list to store all the variant calls. Each variant
    # call is a genomic interval
    variant_interval = []
    for line in infile:
        if (line[0] != '#'):
            pline = line.strip()
            sline = pline.split('\t')
            chrom = 'chr'+sline[7]
            start = int(sline[10])
            end = int(sline[13])

            # Create a 'Genomic interval' from this variant call
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            variant_interval.append(iv)
    infile.close()
    return(variant_interval)





chrom_size_file = open("/home/ubuntu/GRCh38_chrom_size.txt",'r')

# Read chrom size information from the chrom_size_file. 
chrom_size = {}
for line in chrom_size_file:
    pline = line.strip()
    sline = pline.split('\t')
    chrom_size[sline[1]] = int(sline[0])

chrom_size_file.close()

# Creat a 'Genomic Array' using HTSeq package
ga = HTSeq.GenomicArray( chrom_size, stranded=False, typecode="i" )

# Read CNV information from the merged file 
variant_interval = VariantCallTabReader(filepath)

# Get the count of variant calls in each region
variant_num = len(variant_interval)
print "There are "+str(variant_num)+" variant calls from the clustersed studies..."
for i in xrange(variant_num):
    iv = variant_interval[i]
    ga[iv] += 1

outfile = "/home/ubuntu/Yan/test_output.bedgraph"
ga.write_bedgraph_file(outfile, strand=".", track_options="")













