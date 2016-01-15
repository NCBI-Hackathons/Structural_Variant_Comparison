""" Load in the dbVar data and output summary
statistics.

NCBI hackathon dbVar group 2016
"""

import warnings
import pickle, sys
import pandas as pd
import numpy as np
from IPython import embed
import ConfigParser
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from utils import (
        filter_by_size,
        remove_singleton_exp_variants, 
        reverse_dictionary,
        generate_unique_mapping,
        groupby_study_numba,
        )
from generate_report import generate_report

#warnings.filterwarnings("ignore", category=matplotlib.UserWarning)



class FuzzyData(object):
    """
    """
    def __init__(self, data, summary=None):
        self.data = data
        self.summary = summary

    def plot_histogram(self, rpath):
        fig, ax = plt.subplots()
        diff = (self.data.inner_start - 
                self.data.outer_start)
        ax = sns.kdeplot(diff)
        print(diff.mean())
        print(min(diff), max(diff))
        fig.savefig(rpath + 'fuzz_kdeplot.png')







def main(): 
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    report_dict = {}
    gpath = config.get('output', 'output_dir') 
    reader = pd.read_csv(gpath + 'sorted_info.txt', sep="\t", 
            index_col=0, chunksize=500000)
    pool = mp.Pool(4)
    # Begin filtering
    size_limit = config.getfloat('params', 'max_size')
    # Remove duplicated elements
    filtered = []
    func_list = []
    for df in reader:
        f = pool.apply_async(filter_by_size, [df], 
                {'max_size':size_limit})
        func_list.append(f)
    for f in func_list:
        filtered.append(f.get(timeout = 1600))
    df = pd.concat(filtered)
    print(df.shape)
    # :TODO if sstart and sstop are the same, no
    # matter if it was originally annotated as inner_start
    # or inner stop it will be collapsed
    '''
    dfd = df.drop_duplicates(['chr', 'var_type',
        'inner_start', 'start', 'outer_start', 
        'inner_stop', 'stop', 'outer_start'],
        inplace=False)
    '''
    # For now since, for w/e reason
    dfd = df.drop_duplicates(['chr', 'var_type',
        'sstart', 'sstop'])
    new_unique_index = np.arange(dfd.shape[0])
    dfd.loc[:,'uID'] = new_unique_index
    print('new index created')
    # Save intermediate files for now 
    df = generate_unique_mapping(dfd, df, nstudies=2)
    dfd.to_pickle(gpath + 'drop_duplicats_size.pkl')
    df.to_pickle(gpath + 'full_size.pkl')


def study_filtering():
    import timeit
    report_dict = {}
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    gpath = config.get('output', 'output_dir') 
    rpath = config.get('output', 'report_dir')
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
    type_count.round(2)
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

