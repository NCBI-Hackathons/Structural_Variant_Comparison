file1 = '/data/dbVar/Homo_Sapiens/by_study/estd49_Gusev_et_al_2009/tab/estd49_Gusev_et_al_2009.variant_call.submitted.germline.tab.gz'

import gzip
file1 = gzip.open(file1)
next(file1)
line = next(file1).split('\t')
line2= next(file1).split('\t')
for i in range(len(line)):
    print i,line[i],line2[i]
