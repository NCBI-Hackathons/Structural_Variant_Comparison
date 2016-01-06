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
from bx.intervals.intersection import Intersecter, Interval
import seaborn as sns

from utils import copy_test
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
            index_col=0)
    print('Data loaded')
    study_dict = pickle.load(
            open(gpath + 'study_id.dict', 'rb'))
    print('**** study dict loaded ******')
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
    fzd = fuzz_charterize(dfd)
    fz = fuzz_charterize(df)
    fuzz_out = pd.DataFrame({'all_counts': fz.summary['counts'],
        'all_percent': fz.summary['var_percent'],
        'unique_var_counts': fzd.summary['counts'], 
        'unique_var_percent': fzd.summary['var_percent']},
        index = fz.summary.index)
    type_count.to_csv(rpath + 'type_count.txt', sep="\t")
    fuzz_out.to_csv(rpath + 'fuzz_out.txt', sep="\t")
    size_limit = config.get('params', 'size_limit')

    new_unique_index = ['DSV{0!s}'.format(i) for i\
            in xrange(0, dfd.shape[0])]
    dfd['uID'] = new_unique_index
    # Copy number test

    outdf = filter_by_number_of_hits(dfd, df, nstudies=2)
    cnv = copy_test(dfd)
    cnv = cnv.ix[cnv.size <= size_limit, :]


    fig, ax = plt.subplots()
    generate_report(report_dict)
    embed()



if __name__ == '__main__':
    main()
