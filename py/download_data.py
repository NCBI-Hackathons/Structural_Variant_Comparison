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
	f.cwd('pub/dbVar/data/Homo_sapiens/by_study/')
	ldir = list()
	f.retrlines('NLST',ldir.append)
	f.cwd(ldir[0]+'/tab')
	for dir in ldir[1:]:
		ldir2 = list()
		#f.retrlines('NLST',ldir2.append)
		#f.cwd(dir+'/tab')
		f.retrlines('NLST',ldir2.append)
		for dir2 in ldir2:
			#print dir3,type(dir3)
			if 'germline.tab.gz' in dir2 and 'variant_call.' in dir2:
				print 'downloading',dir2
				f2 = open(dir2,'wb')
				f.retrbinary('RETR '+dir2,f2.write)
		#f.cwd('../')
		f.cwd('../../'+dir+'/tab')
	print 'Done'

main()
