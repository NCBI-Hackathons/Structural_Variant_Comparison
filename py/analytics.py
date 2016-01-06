""" Load in the dbVar data and output summary
statistics.

NCBI hackathon dbVar group 2016
"""

import pickle
import pandas as pd
import numpy as np
from IPython import embed
import ConfigParser

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from utils import (filter_by_number_of_hits2,
        remove_singleton_exp_variants, 
        copy_test,
        reverse_dictionary)
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



def collide_intervals(df, fz):
    """
    Arguments
    =========
    df - pandas.dataframe
    fz - fuzzy interval dataframe
    """
    tree_full = Intersecter()
    tree_start = Intersecter()
    fz.apply(lambda x: tree_full.add_interval(x['chr'], 
        x['inner_start'], x['outer_stop']),  
            axis=1)
    fz.apply(lambda x: tree_start.add_interval(x['chr'], 
        x['inner_start'], x['outer_stop']),  
            axis=1)
    return(collision)


def study_chracterize(df):
    """
    """
    pass


def pca_fig(rpath):
    fig, ax = plt.subplots()
    ax.scatter([0,2], [2,6,])
    fig.savefig(rpath + 'test.png')


def fuzz_charterize(df):
    """
    """
    fuzz = np.logical_not(np.logical_or(
        df.outer_start.isnull(),
        df.inner_start.isnull()
        ))
    fz = df.ix[fuzz,:]
    print('Percent fuzzines: {0}'.format(
        float(fuzz.sum())/df.shape[0]))
    by_type = fz.groupby('var_type').agg(lambda x:
        x.shape[0]).loc[:, ['chr']]
    var_percent = by_type.ix[:,0]/float(fz.shape[0])*100
    by_type['var_percent'] = var_percent
    by_type.columns = ['counts', 'var_percent']
    fd = FuzzyData(fz, summary=by_type)
    return(fd)
 

def main(): 
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    report_dict = {}
    gpath = config.get('output', 'output_dir') 
    rpath = config.get('output', 'report_dir')
    df = pd.read_csv(gpath + 'sorted_info.txt', sep="\t", 
            index_col=0, nrows=400000)
    print('Data loaded')
    # Remove duplicated elements
    dfd = df.drop_duplicates(['chr', 
        'inner_start', 'start', 'outer_start', 
        'inner_stop', 'stop', 'outer_start'],
        inplace=False)
    type_count = dfd.groupby('var_type').agg(lambda x:
            x.shape[0]).loc[:, ['chr']]
    var_percent = type_count.ix[:,0]/float(dfd.shape[0])*100
    type_count['var_percent'] = var_percent
    diff = df.shape[0] - dfd.shape[0]
    #tree = Intersecter()
    print('Number of exact duplicate entries: {0}'.
            format(diff))
    print('Unique variant entries: {0}'.format
            (dfd.shape[0]))
    print('Number of records: {0}'.format(df.shape[0]))
    print(type_count)
    report_dict['type_counts'] = type_count.to_html()
    # Get params from file
    size_limit = config.get('params', 'max_size')
    nstudies = config.get('params', 'nstudies')
    new_unique_index = ['DSV{0!s}'.format(i) for i\
            in xrange(0, dfd.shape[0])]
    dfd.loc[:,'uID'] = new_unique_index
    print('new index created')
    # Copy number test
    outdf, udf = filter_by_number_of_hits2(dfd, df, nstudies=2,
            max_size = size_limit)
    #cnv = copy_test(dfd)
    #cnv = cnv.ix[cnv.size <= size_limit, :]
    groups = udf.groupby('var_type')
    from plot import plot_dists
    for name, group in groups:
        plot_dists(group.sstop - group.sstart, name,
                rpath)

    generate_report(report_dict)
    study_dict = pickle.load(
            open(gpath + 'dict_test.txt', 'rb'))
    sdict = reverse_dictionary(study_dict)
    print('**** study dict loaded ******')
    gs = remove_singleton_exp_variants(outdf, sdict,
            nstudies)
    filtered_data = udf.ix[gs.values,:]
    filtered_data.to_csv(gpath + 'filtered_all.txt')


if __name__ == '__main__':
    main()
