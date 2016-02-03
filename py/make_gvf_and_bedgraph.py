#!/bin/bash
# Authors: Yan Kai, Jeffery Hsu
# GRCh38
import sys,re,os
from optparse import OptionParser
import pandas as pd 
import HTSeq

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

def VariantCallTabReader(filepath, chrom_size):
    """
    This function aims to read the variant calls from the merged studies which
    contain various types of variant types. 

    """
    infile = pd.read_csv(filepath, sep="\t")
    
    # var_types is a dic keyed by var_type and valued by a list of genomic intervals
    var_types_ga = {}
    var_types_id = {}
    
    for _, line in infile.iterrows():
        var_type = str(line['var_type'])
        var_type = var_type.replace(" ","_")
        
        if var_type not in var_types_ga.keys():
            var_types_ga[var_type] = []
            var_types_id[var_type] = []
       	
        chrom = 'chr'+ str(line['chr'])
        accesion = line[0]
        if (chrom in chrom_size.keys()):
            start = inner_outer_pref(line, 'start')
            end = inner_outer_pref(line, 'stop')

            # Create a 'Genomic interval' from this variant call
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            var_types_ga[var_type].append(iv)
            var_types_id[var_type].append(accesion)
    return((var_types_ga,var_types_id))

def write_to_gvf(ga,nssd, var_type, outfile):
    """
    ga is the genomic array.
    """
    outfile = open(outfile, 'w')
    ivs = list(ga.steps())
    ids = list(nssd.steps())
    num_iv = len(ivs)
    outline = '##gvf-version 1.00\n'
    outfile.write(outline)
    
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
        SSV = ids[i][1]
        attributes = 'ID='+ID+';'+'count='+str(count)+';'+'SSV='+str(SSV)
        
        outline =chrom+'\t'+source+'\t'+var_type+'\t'+str(start)+'\t'+str(end)+'\t'+score+'\t'+strand+'\t'+phase+'\t'+attributes+'\n'
        outfile.write(outline)
    outfile.close()


def main(argv):
    parser = OptionParser()
    parser.add_option("-r", "--chromsize", action="store", type="string", dest="chromsize", help="GRCh38 chromosome size file", metavar="<str>")
    parser.add_option("-v", "--variantfile", action="store", type="string", dest="variantfile", metavar="<file>", help="the variant calls files in a specific format")
    parser.add_option("-o", "--outdir", action="store", type="string", dest="outdir", metavar="<file>", help="the directory to store the output files")
	
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 6:
        parser.print_help()
        sys.exit(1)

    chrom_size_file = open(opt.chromsize,'r')
    # Read chrom size information from the chrom_size_file. 
    chrom_size = {}
    for line in chrom_size_file:
        pline = line.strip()
        sline = pline.split('\t')
        chrom_size[sline[1]] = int(sline[0])
    chrom_size_file.close()


    var_types_ga = VariantCallTabReader(opt.variantfile,chrom_size)[0]
    var_types_id = VariantCallTabReader(opt.variantfile,chrom_size)[1]
    for var_type in var_types_ga.keys():  
        # Creat a 'Genomic Array' using HTSeq package
        ga = HTSeq.GenomicArray( chrom_size, stranded=False, typecode="i" )
        nssd = HTSeq.GenomicArrayOfSets(chrom_size, stranded=False)
     
        variant_interval = var_types_ga[var_type]
        variant_id = var_types_id[var_type]

        # Get the count of variant calls in each region
        variant_num = len(variant_interval)
        print "For "+var_type+", there are "+str(variant_num)+" variant calls from the clustersed studies..."
        for i in xrange(variant_num):
            iv = variant_interval[i]
            try:
                ga[iv] += 1
                nssd[iv] += variant_id[i]
            except:
                iv.length == 0

        bedgraph = opt.outdir+'/'+var_type+'_dbVar.bedgraph'
        ga.write_bedgraph_file(bedgraph, strand=".", track_options="")
    
        gvf = opt.outdir+'/'+var_type+'_dbVar.gvf'
        write_to_gvf(ga,nssd,var_type,gvf)

if __name__ == "__main__":
    main(sys.argv)
