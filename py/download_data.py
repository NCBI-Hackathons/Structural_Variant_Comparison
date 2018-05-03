#!/usr/bin/env python

#import urllib2
from ftplib import FTP
import os

def main():
        try:
                os.mkdir('dbvar')
        except:
                pass
        os.chdir('dbvar')
        f = FTP('ftp.ncbi.nlm.nih.gov')
        f.login()
        f.cwd('pub/dbVar/data/Homo_sapiens/by_study/tsv/')
        ldir = list()
        f.retrlines('NLST',ldir.append)
        for file in ldir:
                # note: this excludes variant_call.somatic.tsv
                if 'variant_call.tsv.' in file:
                        print 'downloading',file
                        f2 = open(file,'wb')
                        f.retrbinary('RETR '+file,f2.write)
        print 'Done'

main()
