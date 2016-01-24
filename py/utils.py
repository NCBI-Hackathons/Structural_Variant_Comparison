""" Utility functions for filtering
""" 
import numpy as np
import pandas as pd
from bx.intervals.intersection import Intersecter, Interval
import numba


def fuzzy_ends(df):
    sstart = df.loc[:, ['outer_start', 'start','inner_start']].apply(np.min,
            axis=1)
    sstop = df.loc[:, ['outer_stop', 'stop','inner_stop']].apply(np.max, axis=1)
    df.loc[:, 'sstop'] = sstop.astype(np.int32)
    df.loc[:, 'sstart'] =  sstart.astype(np.int32)
    return(df)


def filter_by_size(df, study, max_size=3e6):
    '''
    '''
    df = fuzzy_ends(df)
    df = df.query('sstop - sstart < {0}'.format(float(max_size)))
    s = np.repeat(study , df.shape[0]).astype('|S10')
    df.loc[:, 'study'] = s
    return(df)


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


@numba.jit
def generate_unique_mapping_numba(df_st, df_sp, udf_st,
        udf_sp, dfd_nix):
    """  
    df_st - full dataframe singular start
    df_sp - full dataframe singular stop
    """
    n = len(df_st)
    un = len(dfd_nix)
    out_index = np.zeros(n, dtype=np.int32)
    # Maybe use interval tree here and do it in cython?
    ci = 0
    for i in range(n):
        for j in range(ci, un):
            if (df_st[i] == udf_st[j]) and (df_sp[i] == udf_sp[j]):
                out_index[i] = dfd_nix[j]
                # Just to be safe, assume sorted
                ci = j - 4
                break
            else:pass
    return(out_index)



def generate_unique_mapping(udf, df):
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
            print('KeyError generate_unique_mapping')
            pass
    #df.loc[:, 'uID'] = hits_uids
    return(pd.Series(hits_uids, index=df.index))


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
            print('Yes')
        else:
            more_than_one.append(False)
    out_s = pd.Series(more_than_one, index = uid_index)
    return(out_s, study_list)


def fuzzy_matches(df, stop):
    """
    """
    pass



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
