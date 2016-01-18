#!/usr/bin/env python
""" Create variant_id : study dictionary
"""


import ConfigParser
import pickle
import glob


def get_id(str1):
    name = str1.split('/')[-3]
    return name.split('_')[0]


def main():
    config = ConfigParser.RawConfigParser()
    try:
        config.read('../example.cfg')
        dbVar_path = config.get('input', 'dbVar')
        output_file = config.get('output','output_dir') + 'dict_test.pkl'
    except:
        import sys
        sys.exit('Could not locate config file')
    files = glob.glob(dbVar_path + '/*/tab/*variant_call.*.germline.tab.gz')

    id_dict = {}
    for file_ in files:
        print file_
        study = get_id(file_)
        read = open(file_)
        for line in read:
            if line[0]!='#':
                id_dict[line.split('\t')[0]] = study
                #value.append(line.split('\t')[0])
                #id_dict[key]=value
    save = open(output_file,'w')
    pickle.dump(id_dict,save)
    save.close()


if __name__ == '__main__':
    main()
