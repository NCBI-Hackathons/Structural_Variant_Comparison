output_file = 'dict_test.txt'

import pickle
import glob

def get_id(str1):
    name = str1.split('/')[-3]
    return name.split('_')[0]
def main():
    files=glob.glob('/data/dbVar/Homo_Sapiens/by_study/*/tab/*variant_call.*.germline.tab.gz')
    id_dict = {}
    for file1 in files:
        print file1
        key = get_id(file1)
        read = open(file1)
        value = []
        for Line in read:
            if Line[0]!='#':
                value.append(Line.split('\t')[0])
                id_dict[key]=value
    save = open(output_file,'w')
    pickle.dump(id_dict,save)

main()
