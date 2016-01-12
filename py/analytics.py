""" Load in the dbVar data and output summary
statistics.

NCBI hackathon dbVar group 2016
"""

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
        generate_unique_mapping
        )
from generate_report import generate_report



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




def pca_fig(rpath):
    fig, ax = plt.subplots()
    ax.scatter([0,2], [2,6,])
    fig.savefig(rpath + 'test.png')





def main(): 
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    report_dict = {}
    gpath = config.get('output', 'output_dir') 
    rpath = config.get('output', 'report_dir')
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
    '''
    dfd = df.drop_duplicates(['chr', 'var_type',
        'inner_start', 'start', 'outer_start', 
        'inner_stop', 'stop', 'outer_start'],
        inplace=False)
    '''
    # For now since, for w/e reason
    dfd = df.drop_duplicates(['chr', 'var_type',
        'sstart', 'sstop'])
    new_unique_index = ['DSV{0!s}'.format(i) for i\
            in xrange(0, dfd.shape[0])]
    dfd.loc[:,'uID'] = new_unique_index
    print('new index created')
    # Save intermediate files for now 
    dfd.to_pickle(gpath + 'drop_duplicats_size.pkl')
    df.to_pickle(gpath + 'full_size.pkl')


def study_filtering():
    """ Temp to figure out what is going on
    """
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    gpath = config.get('output', 'output_dir') 
    nstudies = config.get('params', 'nstudies')
    df = pd.read_pickle(gpath + 'full_size.pkl')
    dfd = pd.read_pickle(gpath + 'drop_duplicats_size.pkl')
    df = generate_unique_mapping(dfd, df, nstudies=2)
    print('begin filtering by study')
    study_dict = pickle.load(
            open(gpath + 'dict_test.txt', 'rb'))
    sdict = reverse_dictionary(study_dict)
    print('**** study dict loaded ******')
    gs, sl = remove_singleton_exp_variants(df, sdict,
            nstudies)
    filtered_data = dfd.ix[gs.values,:]
    filtered_data.to_csv(gpath + 'filtered_all.txt')
    embed()


def reports():
    # Begin report generation
    groups = dfd.groupby('var_type')
    from plot import plot_dists
    for name, group in groups:
        plot_dists(group.sstop - group.sstart, name,
                rpath)
    generate_report(report_dict)
    type_count = dfd.groupby('var_type').agg(lambda x:
            x.shape[0]).loc[:, ['chr']]
    var_percent = type_count.ix[:,0]/float(dfd.shape[0])*100
    type_count['var_percent'] = var_percent
    print(type_count)
    report_dict['type_counts'] = type_count.to_html()
    # Get params from file


if __name__ == '__main__':
    if sys.argv[1]=='main':
        main()
    elif sys.argv[1] == 'temp': 
        study_filtering()
    else:
        pass

