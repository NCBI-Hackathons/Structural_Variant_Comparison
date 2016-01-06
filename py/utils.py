""" Pass
""" 
import numpy as np
import pandas as pd
from bx.intervals.intersection import Intersecter, Interval
from collections import defaultdict


def filter_by_number_of_hits(udf, df, 
        ncollisions=3, nstudies=2):
    """ Collide a bx interval tree with a dataframe
    """
    tree = Intersecter()
    start_list = ['outer_start', 'start', 'inner_start']
    stop_list = ['outer_stop', 'stop', 'inner_stop']
    udf['sstart'] = udf.loc[start_list].apply(np.min, axis=1)
    udf['sstop'] = udf.loc[stop_list].apply(np.max, axis=1) 
    # :TODO  parallelize this
    tree_dict = {}
    for i, j in udf.iterrows():
        tree_dict[j['chr']].add_interval(Interval(j.sstart, j.send, 
            value={'id':j.uID}))

    def _apply(x, ):
        """Stuff
        """

    for i, j in df.iterrows():
        hits = tree_dict[j['chr']].find()
        for 


def copy_test(df, stuy_dict=None):
    """ Copy number variant testing
    """
    # :TODO change to groups
    dfg = ((df['var_type'] == 'copy number gain') |\
           ( df['var_type'] == 'copy number loss') |\
           ( df['var_type'] == 'copy number variation'))
    # append together
    cnv = df.ix[dfg, :]
    def _get_sizes(df):
        start = df.ix[:,['outer_start', 'start',
            'inner_start']].apply(np.min, axis=1)
        end = df.ix[:, ['inner_stop', 'stop',
            'outer_stop']].apply(np.min, axis=1)
        diff = end - start
        return(diff)
    diff_cnv = _get_sizes(cnv)
    cnv['size'] = diff_cnv
    return(cnv) 
