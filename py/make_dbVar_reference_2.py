import numpy as np
import pandas as pd 
import math
import HTSeq
#filepath = "/home/ubuntu/Yan/estd211_Campbell_et_al_2014.variant_call.remap.germline.tab"
filepath = '/home/ubuntu/dbvar_data/sorted_info.txt'
output_dir = '/home/ubuntu/dbvar_data/output'

def inner_outer_pref(line, suffix):
    ''' 
    Preferences outer + suffix -> suffix -> inner + suffix
    '''
    prefix = ['inner_', '', 'outer_']
    test = [i + suffix for i in prefix]
    for i in test:
        if(not pd.isnull(line[i])):
            out = int(line[i])
    return(out)


def VariantCallTabReader(filepath,chrom_size):
    """
    1. This function is to read from the merged file to obatin the genomic
    regions of each variant call. 
    2. Variant type is assumed to be the same (i.e only deletion at this
    stage). 
    3. No fussiness infomation are considered.

    """
    infile = pd.read_csv(filepath, sep="\t")

    # variant_interval is a list to store all the variant calls. Each variant
    # call is a genomic intervaln
    variant_interval = []
    for _, line in infile.iterrows():
        chrom = 'chr'+ str(line['chr'])
        if (chrom in chrom_size.keys()):
            start = inner_outer_pref(line, 'start')
            end = inner_outer_pref(line, 'stop')

         # Create a 'Genomic interval' from this variant call
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            variant_interval.append(iv)
    return(variant_interval)
    
def VariantCallTabReader_V2(filepath, chrom_size):
    """
    This function aims to read the variant calls from the merged studies which
    contain various types of variant types. 

    """
    infile = pd.read_csv(filepath, sep="\t")
    
    # var_types is a dic keyed by var_type and valued by a list of genomic intervals
    var_types_ga = {}
	
    for _, line in infile.iterrows():
        
        var_type = str(line['var_type'])
        
        if var_type not in var_types_ga.keys():
            var_types_ga[var_type] = []
           	
        	
        chrom = 'chr'+ str(line['chr'])
        if (chrom in chrom_size.keys()):
            start = inner_outer_pref(line, 'start')
            end = inner_outer_pref(line, 'stop')

            # Create a 'Genomic interval' from this variant call
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            var_types_ga[var_type].append(iv)
    return(var_types_ga)


def write_to_gvf(ga,var_type, outfile):
    """
    ga is the genomic array.
    """
    outfile = open(outfile, 'w')
    ivs = list(ga.steps())
    num_iv = len(ivs)
    outline = '##gvf-version 1.00'

    source = 'dbVar'
    score = '.'
    strand = '.'
    phase = '.'
    for i in xrange(num_iv):
        chrom = ivs[i][0].chrom
        start = ivs[i][0].start+1
        end = ivs[i][0].end
        ID = str(i+1)

        count = ivs[i][1]
        attributes = 'ID='+ID+';'+'count='+str(count)
        outline =chrom+'\t'+source+'\t'+var_type+'\t'+str(start)+'\t'+str(end)+'\t'+score+'\t'+strand+'\t'+phase+'\t'+attributes+'\n'
        outfile.write(outline)
    outfile.close()




chrom_size_file = open("/home/ubuntu/GRCh38_chrom_size.txt",'r')

# Read chrom size information from the chrom_size_file. 
chrom_size = {}
for line in chrom_size_file:
    pline = line.strip()
    sline = pline.split('\t')
    chrom_size[sline[1]] = int(sline[0])

chrom_size_file.close()



var_types_ga = VariantCallTabReader_V2(filepath,chrom_size)
for var_type in var_types_ga.keys():
    
    # Creat a 'Genomic Array' using HTSeq package
    ga = HTSeq.GenomicArray( chrom_size, stranded=False, typecode="i" )
     
    variant_interval = var_types_ga[var_type]

    # Get the count of variant calls in each region
    variant_num = len(variant_interval)
    print "For "+var_type+", there are "+str(variant_num)+" variant calls from the clustersed studies..."
    for i in xrange(variant_num):
        iv = variant_interval[i]
        try:
            ga[iv] += 1
        except:
            iv.length == 0

    bedgraph = output_dir+'/'+var_type+'_dbVar.bedgraph'
    ga.write_bedgraph_file(bedgraph, strand=".", track_options="")

    gvf = output_dir+'/VarType='+var_type+'_dbVar.gvf'
    write_to_gvf(ga,var_type,gvf)


