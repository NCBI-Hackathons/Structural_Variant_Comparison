""" 
TriLe
"""

import gzip
import glob
import pickle

def ReadTab(input_file,output_file):
    file1 = gzip.open(input_file)
    save = open(output_file,'a')
    next(file1)
    next(file1)
    for col in file1:
        col = col.split('\t')
        if col[6]=='GRCh38':
            save.write('\t'.join([col[0],col[2],col[7],col[9],col[10],col[11],col[12],col[13],col[14],col[45]])+'\n')
    save.close()


def main(output_file):
    files = glob.glob('/data/dbVar/Homo_Sapiens/by_study/*/tab/*variant_call.*.germline.tab.gz')
    save = open(output_file,'w')
    save.write('#accession_num\tvar_type\tchr\touter_start\tstart\tinner_start\tinner_stop\tstop\touter_stop\tremap_alignment\n')
    save.close()
    for i in files:
        print 'Reading',i
        ReadTab(i,output_file)

if __name__ == '__main__':
    main('unsorted_info.txt')
