#!/usr/bin/env python
""" Load in the dbVar data and output summary
statistics.

NCBI hackathon dbVar group 2016
"""

import timeit, glob
import pickle, sys
import pandas as pd
import numpy as np
import ConfigParser
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from utils import (
        filter_by_size,
        reverse_dictionary,
        generate_unique_mapping,
        groupby_study_numba,
        generate_unique_mapping_numba,
        )
from generate_report import generate_report

from IPython import embed


def append_p(study, x):
    study = study.split("/")[-1].rstrip(".txt")
    s=np.array([study for i in range(x.shape[0])], dtype='|S10')
    x.loc[:,'study'] = s
    return(x)



def main(): 
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    gpath = config.get('output', 'output_dir') 
    size_limit = config.getfloat('params', 'max_size')
    files = glob.glob("/home/ubuntu/gvf_by_studies/tab/*.txt")
    studies_include = config.get('params', 'studies_include')
    studies_exclude = config.get('params', 'studies_exclude').split(",")
    if studies_include == '' or studies_include == None:
        studies_include = []
    else:
        studies_include = studies_include.split(",")
    filtered = []
    start = timeit.default_timer()
    pool = mp.Pool(8)
    files = files[0:10]
    studies = [i.split("/")[-1].rstrip(".txt") for i in files]
    for i in files:
        study = i.split("/")[-1].rstrip(".txt")
        if study in studies_exclude: pass
        else:
            if (len(studies_include) == 0) or (study in studies_include):
                print('yeah')
                reader = pd.read_csv(i, sep="\t", 
                        index_col=0, dtype={'chr':'S5'})
                pool.apply_async(filter_by_size, [reader, study],
                        {'max_size': size_limit},
                        callback = lambda x: filtered.append(x))
            else: pass
    # Remove duplicated elements
    ###### Step takes around 7 minutes ###################
    pool.close()
    pool.join()
    df = pd.concat(filtered)
    stop = timeit.default_timer()
    print('Time to load in files and parse: {0!s}'.format(stop-start))
    p_studies = set(df.study)
    non_passed = []
    for i in studies:
        if i not in p_studies:
            non_passed.append(i)
    print('Studies that had no variants that did not pass size filtering:{0}'.format("\t".join(non_passed)))
    ############## HACK for now until we find out what is going on #
    df = df.ix[np.logical_not(df.index.duplicated()),:]
    # :TODO if sstart and sstop are the same, no
    # matter if it was originally annotated as inner_start
    # or inner stop it will be collapsed
    # For now since, ignore fuzzy 
    dfd = df.drop_duplicates(['chr', 'var_type',
        'sstart', 'sstop'], inplace=False)
    new_unique_index = np.arange(dfd.shape[0])
    dfd.loc[:,'uID'] = new_unique_index
    print('new index created')
    # This step takes forever
    start = timeit.default_timer()
    groups = df.groupby('chr')
    unique_mapping = []
    pool = mp.Pool(8)
    for name, group in groups:
        pool.apply_async(generate_unique_mapping,
                args = (dfd.ix[dfd.chr == name,:], group),  
                callback=lambda x: unique_mapping.append(x))
        '''
        tgroup = dfd.ix[dfd['chr'] == name,]
        pool.apply_async(generate_unique_mapping_numba,
                args = (group.sstart.values, 
                    group.sstop.values, 
                    tgroup.sstart.values, 
                    tgroup.sstop.values, 
                    tgroup.index.values),
                callback=lambda x: unique_mapping.append(pd.Series(x,
                    index = group.index)))
        '''
    pool.close()
    pool.join()
    ns = pd.concat(unique_mapping)
    stop = timeit.default_timer()
    embed()
    print('Time to generate mapping: {0!s}'.format((stop-start)))
    df['uID'] = ns


def study_filtering():
    report_dict = {}
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    gpath = config.get('output', 'output_dir') 
    rpath = config.get('output', 'report_dir')
    vartype_f = config.('params', 'var_type')
    nstudies = config.getint('params', 'nstudies')
    df = pd.read_pickle(gpath + 'full_size.pkl')
    dfd = pd.read_pickle(gpath + 'drop_duplicats_size.pkl')
    print('begin filtering by study')
    study_dict = pickle.load(
            open(gpath + 'dict_test.txt', 'rb'))
    sdict = reverse_dictionary(study_dict)
    print('**** study dict loaded ******')
    s1 = np.array([sdict[i] for i in df.index], dtype='|S20')
    start = timeit.default_timer()
    output = np.zeros(dfd.uID.shape[0], dtype=bool)
    std_filter = groupby_study_numba(df.uID.values, s1, 
            output, nstudies=nstudies) 
    stop = timeit.default_timer()
    print(np.sum(std_filter))
    dfd = dfd.ix[std_filter,:]
    df = df.ix[df.uID.isin(dfd.uID),:]
    dfd.to_csv(gpath + 'filtered_no_dupes.txt', sep="\t")
    df.to_csv(gpath + 'study_filtered_all.txt', sep="\t")
    print('Time to run: {0!s}'.format(stop - start))
    groups = dfd.groupby('var_type')
    from plot import plot_dists
    generate_report(report_dict)
    for name, group in groups:
        plot_dists(group.sstop - group.sstart, name,
                rpath)
    type_count = dfd.groupby('var_type').agg(lambda x:
            x.shape[0]).loc[:, ['chr']]
    var_percent = type_count.ix[:,0]/float(dfd.shape[0])*100
    type_count['var_percent'] = var_percent
    type_count['var_percent'].round(2)
    report_dict['var_type_pivot'] = type_count.to_html()
    report_dict['studies'] = []
    report_dict['var_types'] = [name for name, _ in groups]
    generate_report(report_dict)





if __name__ == '__main__':
    if sys.argv[1]=='main':
        main()
    elif sys.argv[1] == 'temp': 
        study_filtering()
    else:
        pass

