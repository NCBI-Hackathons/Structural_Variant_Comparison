""" Pass
""" 
import numpy as np
import pandas as pd
from bx.intervals.intersection import Intersecter, Interval
from collections import defaultdict
from IPython import embed


def fuzzy_ends(df):
    tdf = df.loc[:, ['outer_start', 'start','inner_start']]
    sstart = tdf.apply(np.min, axis=1)
    df.loc[:,'sstart'] =  sstart
    tdf = df.loc[:, ['outer_stop', 'stop','inner_stop']]
    sstop = tdf.apply(np.max, axis=1)

    df.loc[:, 'sstop'] = sstop
    return(df)




    

def fuzzy_diff(df):
    """ This was the third try for this, still
    very slow. 
    """
    sstart = np.zeros(df.shape[0], dtype=np.int64)
    sstop = np.zeros(df.shape[0], dtype=np.int64)
    enum = 0
    for _, j in df.iterrows():
        if not pd.isnull(j['outer_start']):
            sstart[enum] = j['outer_start']
        else:
            try:
                sstart[enum] = j['start']
            except ValueError:
                sstart[enum] = j['inner_start']
        if not pd.isnull(j['outer_stop']):
            sstop[enum] = j['outer_stop']
        else:
            try:
                sstop[enum] = j['stop']
            except ValueError:
                sstop[enum] = j['inner_stop']
        enum += 1
    df.loc[:,'sstart'] =  sstart
    df.loc[:,'sstop'] = sstop
    return(df)


def filter_by_number_of_hits2(udf, df, 
        ncollisions=3, nstudies=2, max_size=3e6):
    ''' Exact matching by sstop and sstart
    '''
    print('starting size filtering')
    udf = fuzzy_ends(udf)
    diff = (udf.sstop - udf.sstart)
    udf = udf.ix[diff < float(max_size), :]
    print('flitered unique')
    df = fuzzy_ends(df)
    diff_full = (df.sstop - df.sstart)
    df = df.ix[diff_full < float(max_size),:]
    print('finished filtering')
    hits_uids = []
    comp_dict = {}
    for _, j in udf.iterrows():
        comp_dict[(j['chr'], j.sstart, j.sstop)] = j.uID
    print('finished generating uID dictionary')
    for _, j in df.iterrows():
        hits_uids.append(comp_dict[(j['chr'],j.sstart, j.sstop)])
        '''
        matches = ((udf.sstart == j.sstart) &
                (udf.sstop == j.sstop))
        match_ids = udf.ix[matches, 'uID']
        if len(match_ids):
            hits_uids.append(list(match_ids.values)[0])
        else:
            hits_uids.append(list(match_ids.values)[0])
        '''
    df['uID'] = hits_uids
    return(df, udf)


def reverse_dictionary(dictionary):
    new_dict = {}
    for i, j in dictionary.iteritems():
        for k in j:
            new_dict[k] = i
    return(new_dict)
        

def remove_singleton_exp_variants(df, study_dict,
        nstudies=2):
    """
    """
    ugroups = df.groupby('uID')
    more_than_one = []
    uid_index = []
    for name, group in ugroups:
        studies = [study_dict[i] for i in group.index]
        studies = set(studies)
        uid_index.append(name)
        if len(studies) > nstudies:
            more_than_one.append(True)
        else:
            more_than_one.append(False)
    out_s = pd.Index(more_than_one, index = uid_index)
    return(out_s)



def copy_test(df):
    """ Copy number variant testing
    """
    # :TODO change to groups
    dfg = ((df['var_type'] == 'copy number gain') |\
           ( df['var_type'] == 'copy number loss') |\
           ( df['var_type'] == 'copy number variation'))
    # append together
    cnv = df.ix[dfg, :]
    #cnv['size'] = cnv.sstop - cnv.sstart
    return(cnv) 
