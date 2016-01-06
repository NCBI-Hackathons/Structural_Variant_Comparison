#!/usr/bin/env python
#
#
USAGE = '''
USAGE (all parameters are OPTIONAL):
    -size\ta:b / a: / :b\tRange of variant size
    -chr\tn1 n2 ...\tTargeted chromosome
    -count\ta:b / a: / :b\tRange of variant count
    -range\ta:b / a: / :b\tRange of searching position
    -type\tvar_type1 var_type2 ...\tFor a list variant type, please find it in our manual
    -output\toutput_file\tDefault is "filtered_output.gvf")
\n\n
'''

####
from sys import argv
import os
import ConfigParser
import glob
####


#### Getting database
database = glob.glob('/home/ubuntu/dbvar_data/output/*.gvf')
#print database
####


####
def check_chr(base,var):
    if base == []:
        return True
    var = var[3:]
    if var.lower() in base or var.upper() in base:
        return True
    else:
        return False

def check_count(base,var):
    var = int(var.split(';')[1].split('=')[-1])
    base = base.split(':')
    if base == ['']:
        return True

    if base[0]=='':
        base[0]=0
    if base[1]=='':
        base[1]='inf'

    base[0]=float(base[0])
    base[1]=float(base[1])
    #print base,var

    if var>base[1] or var<base[0]:
        return False
    else:
        return True

def check_type(base,var):
    if base == []:
        return True

    for item in base:
        if item.lower()==var.lower():
            return True
    return False
            
def check_region(base_region,var_region):
    base_region = base_region.split(':')
    if base_region == ['']:
        return True

    if base_region[0]=='':
        base_region[0]=0
    if base_region[1]=='':
        base_region[1]='inf'
    base_start,var_start=float(base_region[0]),float(var_region[0])
    base_end,var_end=float(base_region[1]),float(var_region[1])
    if var_start>base_region:
        return 'FalseWithBreak'
    elif var_start<base_start or var_end>base_end:
        return False
    else:
        return True

def get_data(key,data_type,default=None):
    #Get data from command line
    if key not in argv:
        if default == None:
            return data_type()
        else:
            return default

    key_in = argv.index(key)
    if data_type==str:
        output = argv[key_in+1]
    elif data_type==list:
        output = []
        for i in range(key_in+1,len(argv)):
            if argv[i][0]=='-':
                break
            else:
                output.append(argv[i])
    return output


def check_size(base,var):
    base = base.split(':')
    if base==['']:
        return True

    if base[0]=='':
        base[0]=0
    if base[1]=='':
        base[0]='inf'
    base[0]=float(base[0])
    base[1]=float(base[1])
    var = int(var[1])-int(var[0])+1

    if var>base[1] or var<base[0]:
        return False
    else:
        return True

def main():
    config = ConfigParser.RawConfigParser()
    if '-h' in argv:
        print USAGE
        exit()
    output_file = get_data('-output',str,'filtered_output.gvf')
    base_region = get_data('-range',str)
    var_type = get_data('-type',list)
    base_count = get_data('-count',str)
    base_chr = get_data('-chr',list)
    base_size = get_data('-size',str)
    if output_file == None:
        output_file = os.path.join( 
                config.get('output', 'output_dir'), 
                'filtered_output.gvf')


    save = open(output_file,'w')
    for data in database:
        read = open(data)
        for Line in read:
            if Line[0]!='#':
                line = Line.split('\t')
                var_region = line[3:5]
                key1 = check_region(base_region,var_region)
                key2 = check_type(var_type,line[2])
                key3 = check_count(base_count,line[8])
                key4 = check_chr(base_chr,line[0])
                key5 = check_size(base_size,var_region)
                keys = [key1,key2,key3,key4,key5]

                if 'FalseWithBreak' in keys:
                    break
                elif False not in keys:
                    save.write(Line)
    save.close()
####

main()
