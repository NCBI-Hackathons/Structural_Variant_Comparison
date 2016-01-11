
'''
fz = fuzz_charterize(df)
fuzz_out = pd.DataFrame({'all_counts': fz.summary['counts'],
    'all_percent': fz.summary['var_percent'],
    'unique_var_counts': fzd.summary['counts'], 
    'unique_var_percent': fzd.summary['var_percent']},
    index = fz.summary.index)
type_count.to_csv(rpath + 'type_count.txt', sep="\t")
fuzz_out.to_csv(rpath + 'fuzz_out.txt', sep="\t")
'''

def fuzzy_mapping(udf, df, 
        ncollisions=3, nstudies=2, max_size=3e6):
    """ Collide a bx interval tree with a dataframe

    Fuzzy matching
    """
    start_list = ['outer_start', 'start', 'inner_start']
    stop_list = ['outer_stop', 'stop', 'inner_stop']
    sstart = udf.loc[:,start_list].apply(np.min, axis=1)
    udf.loc[:, 'sstart'] =  sstart
    sstop = udf.loc[:,stop_list].apply(np.max, axis=1)
    udf.loc[:, 'sstop'] = sstop
    # :TODO  parallelize this
    diff = (udf.sstop - udf.start)
    udf = udf.ix[diff < max_size, :]
    tree_dict = {}
    for i in range(1, 23):
        tree_dict[str(i)] = Intersecter()
    tree_dict['X'] = Intersecter()
    tree_dict['Y'] = Intersecter()
    for i, j in udf.iterrows():
        try:
            tree_dict[str(j['chr'])].add_interval(
                    Interval(j.sstart, j.sstop, value={'id':j.uID}
                        ))
        except KeyError:
            pass
    print("######### Finished creating tree ##########")
    #:TODO include all samples
    #:TODO fuzz intervals need to be their own interval trees
    df = df.iloc[0:20,:]
    def _apply_find_exact(x, tree_dict):
        """Stuff
        """
        hits = tree_dict[x['chr']].find(x.sstart, x.sstop)
        if hits >= 1:
            for i in hits:
                if (i.start == x.sstart &\
                        i.end == x.end):
                    return(i.uID)
                else: 
                    return(np.nan)
        # Return only exact matches
    # this is ridiculously low 
    z = df.apply(_apply_find_exact, args=(tree_dict,))
    df['uID'] = z
    return(df)


########### Report generation

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
