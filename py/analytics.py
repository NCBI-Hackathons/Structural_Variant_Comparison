""" Load in the dbVar data and output summary
statistics.

NCBI hackathon dbVar group 2016
"""

import pickle
import pandas as pd
import numpy as np
from IPython import embed
import HTSeq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bx.intervals.intersection import Intersecter, Interval



def collide_intervals(dataframe):
    """
    Arguments
    =========
    dataframe - pandas.dataframe
    """
    raise NotImplemented
    return(collision)


def pca_fig():
    fig, ax = plt.subplots()
    ax.scatter([0,2], [2,6,])
    plt.show()


def main(): 
    gpath = '/home/ubuntu/dbvar_data/sorted_info.txt'
    df = pd.read_csv(gpath, sep="\t")
    print('Data loaded')

    # Remove duplicated elements
    dfd = df.drop_duplicates(['chr', 
        'inner_start', 'start', 'outer_start', 
        'inner_stop', 'stop', 'outer_start'],
        inplace=False)
    type_count = dfd.groupby('var_type').agg(lambda x:
    x.shape[0]).loc[['var_type']]
    var_percent = type_count/float(dfd.shape[0])*100
    type_count['var_percent'] = var_percent
    diff = df.shape[0] - dfd.shape[0]
    print('Number of exact duplicate entries: {0}'.
            format(diff))
    print('Unique variant entries: {0}'.format
            (dfd.shape[0]))
    '''
    bychrom = df.groupby('chr')

    for name, group in bychrom:
        group = group.ix[]
    '''
    fuzz3 = np.logical_not(np.logical_or(df.outer_start.isnull(),
        df.inner_start.isnull()))
    #tree = Intersecter()
    print('Number of records: {0}'.format(df.shape[0]))
    print('Percent fuzzines: {0}'.format(
        float(fuzz3.sum())/df.shape[0]))
    pca_fig()
    embed()



if __name__ == '__main__':
    main()
