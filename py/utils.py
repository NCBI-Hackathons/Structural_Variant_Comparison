""" Pass
""" 
import numpy as np
import pandas as pd
#from bx.intervals.intersection import Intersecter, Interval
from IPython import embed
import numba



def fuzzy_ends(df):
    tdf = df.loc[:, ['outer_start', 'start','inner_start']]
    sstart = tdf.apply(np.min, axis=1)
    df.loc[:,'sstart'] =  sstart.astype(np.int32)
    tdf = df.loc[:, ['outer_stop', 'stop','inner_stop']]
    sstop = tdf.apply(np.max, axis=1)
    df.loc[:, 'sstop'] = sstop.astype(np.int32)
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


def filter_by_size(df, max_size=3e6):
    '''
    '''
    df = fuzzy_ends(df)
    diff = (df.sstop - df.sstart)
    df = df.ix[diff < float(max_size), :]
    return(df)


def generate_unique_mapping(udf, df, 
        ncollisions=3, nstudies=2):
    ''' Exact matching by sstop and sstart
    '''
    hits_uids = []
    comp_dict = {}
    for _, j in udf.iterrows():
        comp_dict[(j['chr'], j['var_type'],
            j.sstart, j.sstop)] = j.uID

    for _, j in df.iterrows():
        try:
            hits_uids.append(comp_dict[(j['chr'],
                 j['var_type'],
                j.sstart, j.sstop)])
        except KeyError:
            # there shouldn't be any key errors
            '''
            qstring = ('sstart >= {0} & '
                    'sstop <= {1} & '
                    'chr == {2}'
                    'var_type == {3}'
                    )
            qtest = udf.query(qstring.format(j.sstart, 
                j.sstop, j.chr))
            embed()
            '''
            print('KeyError generate_unique_mapping')
            pass
    df['uID'] = hits_uids
    return(df)


def reverse_dictionary(dictionary):
    new_dict = {}
    for i, j in dictionary.iteritems():
        for k in j:
            new_dict[k] = i
    return(new_dict)


@numba.jit
def groupby_study_numba(index, value, output, 
        nstudies = 2, gold_standard = None):
    """
    gold_standard - alist of studies with variants that 
    will be ignored when checking for singletons
    """
    # Need to avoid dictionaries for numba 
    sl = np.zeros(len(output), dtype='|S200')
    for i in range(index.shape[0]):
        sl[index[i]] += value[i] + ','
    z = np.char.count(sl, ',')
    return(z >= nstudies)
    




        
        

def remove_singleton_exp_variants(df, study_dict,
        nstudies=2):
    """
    """
    more_than_one = []
    uid_index = []
    study_list = []
    ugroups = df.groupby('uID')
    # Not sure why this isn't working
    for name, group in ugroups:
        studies = [study_dict[i] for i in group.index]
        study_list.extend(studies)
        studies = set(studies)
        uid_index.append(name)
        if len(studies) >= nstudies:
            more_than_one.append(True)
        else:
            more_than_one.append(False)
    out_s = pd.Series(more_than_one, index = uid_index)
    return(out_s, study_list)



def copy_test(df):
    """ Copy number variant testing
    """
    # :TODO change to groups
    dfg = ((df['var_type'] == 'copy number gain') |\
           ( df['var_type'] == 'copy number loss') |\
           ( df['var_type'] == 'copy number variation'))
    # append together
    cnv = df.ix[dfg, :]
    return(cnv) 
