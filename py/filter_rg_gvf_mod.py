#!/usr/bin/env python
import os
import sys
import ConfigParser
import time
import glob

#database = ['VarType=complex_dbVar.gvf']
#database = glob.glob('/home/ubuntu/gvf_by_studies/by_studies/*/*.gvf')
database = glob.glob('/home/ubuntu/output_gvf/*.gvf')
config_header_filter = 'filter'

def get_parameter(key_word,data_type,default):
	cmd = sys.argv
	output = data_type()
	if '-config' in cmd:
		config_file = cmd[cmd.index('-config')+1]
		config = ConfigParser.ConfigParser()
		config.read(config_file)
		try:
			output = config.get(config_header_filter,key_word)
			if data_type == list:
                                output = output.split()
                        else:
                                output = data_type(output)
		except ConfigParser.NoSectionError:
			print '\nContiguration file is in wrong format\n'
			exit()
		except ConfigParser.NoOptionError:
			output = data_type()
	elif '-'+key_word in cmd:
		output = data_type()
		for i in range(cmd.index('-'+key_word)+1,len(cmd)):
			if cmd[i][0]=='-':
				break
			elif data_type == list:
				output.append(cmd[i])
			else:
				output += data_type(cmd[i])

	if output == data_type():
		print 'Set ['+key_word+'] as default'
		return default
	else:
		print 'Set ['+key_word+'] as '+str(output)
		return output

def get_range(string):
	string = string.split(':')
	if string[0]=='':
		start = 0
	else:
		start = int(string[0])
	if string[1]=='':
		end = float('inf')
	else:
		end = int(string[1])
	return [start,end]

def standardize_list(a_list,data_type=str):
	for i in range(len(a_list)):
		if data_type==str:
			a_list[i]==a_list[i].lower()
		else:
			a_list[i]==data_type(a_list[i])
	return a_list

def check_range(search_range,line):
	line = line.split()
	if int(line[3])>search_range[1]:
		return 'FalseWithBreak'
	elif int(line[3])<search_range[0] or int(line[4])>search_range[1]:
		return False
	else:
		return True

def check_count(search_count,line):
	line = line.split()[8].split(';')
	for item in line:
		if 'count' in item:
			count = int(item.replace('\n','').split('=')[-1])
			break
	if count<search_count[0] or count>search_count[1]:
		return False
	else:
		return True

def check_type(search_type,line):
	var_type = line.split()[2]
        if search_type == 'all':
            return True
	if var_type in search_type:
		return True
	else:
		return False

def check_chr(search_chr,line):
	var_chr = line.split()[0][3:]
	if search_chr == ['all']:
		return True
	elif var_chr in search_chr:
		return True
	else:
		return False

def check_len(search_len,line):
	line = line.split()
	var_len = int(line[4])-int(line[3])+1
	if var_len<search_len[0] or var_len>search_len[1]:
		return False
	else:
		return True

def standardize_dir(a_dir):
	if a_dir[-1]=='/':
		a_dir = a_dir[:-1]
	return a_dir

def main():
	print
	#Get parameters
        start = time.time()
	search_range = get_range(get_parameter('range',str,':'))
	search_count = get_range(get_parameter('count',str,':'))
	search_chr = standardize_list(get_parameter('chr',list,['all']))
	search_type = standardize_list(get_parameter('type',list,'all'))
	search_len = get_range(get_parameter('size',str,':'))
	output_dir = standardize_dir(get_parameter('output',str,os.getcwd()))
	#Initilize
	try:
		os.makedirs(output_dir)
		print 'Create new directory',output_dir
	except OSError:
		print 'Overwrite files in directory:',output_dir
        filtered_data = open('/'.join([output_dir,'filtered_data.gvf']),'w')
	#Filter
        print '\nFiltering...'
	for data in database:
		data = open(data)
		for line in data:
			if line[0]!='#':
				status = []
				status.append(check_range(search_range,line))
				status.append(check_type(search_type,line))
				status.append(check_count(search_count,line))
				status.append(check_chr(search_chr,line))
				status.append(check_len(search_len,line))
				if 'FalseWithBreak' in status:
					break
				elif False not in status:
					filtered_data.write(line)
        print 'Done! Running time: '+str(round(time.time()-start,2))+'s\n'
if __name__ == '__main__':
	main()
