#/bin/bash/
'''
This script is for collecting and merging variant based on a given study id
Usage:
    -id [id1 id2 ...] study id
    -output (default is merged_data.txt)
'''

import gzip
import glob
import pickle
from sys import argv

database = glob.glob('/data/dbVar/Homo_Sapiens/by_study/*/tab/*variant_call.*.germline.tab.gz')

def get_data(key,data_type,default=None):
    if key not in argv:
        if default == None:
            return data_type()
        else:
            return default

    key_in = argv.index(key)
    if data_type==str:
        output = argv[key_in+1]
    elif data_type==list:
        output = []
        for i in range(key_in+1,len(argv)):
            if argv[i][0]=='-':
                break
            else:
                output.append(argv[i])
    return output

def read_file(file1,save):
    file1 = gzip.open(file1)
    for Line in file1:
        if Line[0]!='#':
            col = Line.split('\t')
            if col[6]=='GRCh38':
                save.write('\t'.join([col[0],col[2],col[7],col[9],col[10],col[11],col[12],col[13],col[14],col[45]])+'\n')

def main():
    output_file = get_data('-output',str,'merged_data.txt')
    
    all_study_id = []
    for name in database:
        name = name.split('/')[-1].split('_')[0]
        all_study_id.append(name)

    input_study_id = get_data('-id',list,all_study_id)
    
    save = open(output_file,'w')
    for study_id in input_study_id:
        file_in_study = glob.glob('/data/dbVar/Homo_Sapiens/by_study/'+study_id+'*/tab/*variant_call.*.germline.tab.gz')
        for file1 in file_in_study:
            read_file(file1,save)
    save.close()

main()
