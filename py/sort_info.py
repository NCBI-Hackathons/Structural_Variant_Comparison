""" Sort the file created by **
"""
import pandas as pd
from IPython import embed


def TraceBack(key,input_file,output_file):
    save = open(output_file,'a')
    file1 = open(input_file)
    next(file1)
    for Line in file1:
        line = Line.split('\t')
        #key_list = (Int(line[2]),Int(line[3]),Int(line[4]))
        if key==line[0]:
            save.write(Line)
            break
    save.close()


def main(input_file,output_file):
    save = open(output_file,'w')
    save.write('#accession_num\tvar_type\tchr\touter_start\tstart\tinner_start\tinner_stop\tstop\touter_stop\tremap_alignment\n')
    save.close()
    df = pd.read_csv(input_file, sep="\t")
    # :TODO add back chrom
    # outer_start has priority, followed by start then inner_start
    dummy =  df.outer_start.copy(deep=True)
    dummy[pd.isnull(dummy)] = df.start[pd.isnull(dummy)]
    dummy[pd.isnull(dummy)] = df.inner_start[pd.isnull(dummy)]
    assert pd.isnull(dummy).sum() == 0
    df['dummy'] = dummy
    print('sorting')
    df = df.sort_values(by = ['chr', 'dummy'])
    del df['dummy']
    df.to_csv(output_file, sep="\t", index=False)
    '''
    file1 = open(input_file)
    next(file1)
    dict1 = {}
    print 'Sorting...'
    for line in file1:
        line = line.split('\t')
        i = 2
        while line[i]=='':
            i+=1
        start = int(line[i])
        dict1[start]=line[0]
    list1 = dict1.keys()
    list1.sort()
    print 'Tracing back...'
    for item in list1:
        TraceBack(dict1.get(item),input_file,output_file)
    '''

if __name__ == '__main__':
    main('unsorted_info.txt','sorted_info.txt')
