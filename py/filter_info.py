'''
TriLe
'''

USAGE = '''
USAGE:
    -range start end\t(Optional)
    -type variant_type\t(Optional)
    -output output_file\t(Optional - Default is "filtered_output.gvf")
EXAMPLE CMD:
    filter_info.py -range 100 1000 -type CNV Indel dup -output
    my_project/filtered_data2.gvf
'''

####
from sys import argv
import glob
####


#### Getting database
database = glob.glob('/home/ubuntu/Structural_Variant_Comparison/VarType=*.gvf')
#print database
####


####
def check_type(base,var):
    if base == []:
        return True

    for item in base:
        if item.lower()==var.lower():
            return True
    return False
            
def check_region(base_region,var_region):
    if base_region == []:
        return True

    base_start,var_start=int(base_region[0]),int(var_region[0])
    base_end,var_end=int(base_region[1]),int(var_region[1])
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

def main():
    output_file = get_data('-output',str,'filtered_output.gvf')
    base_region = get_data('-range',list)
    var_type = get_data('-type',list)

    save = open(output_file,'w')
    for data in database:
        read = open(data)
        for Line in read:
            line = Line.split('\t')
            var_region = line[3:5]
            key1 = check_region(base_region,var_region)
            key2 = check_type(var_type,line[2])
            keys = [key1,key2]

            if 'FalseWithBreak' in keys:
                break
            elif False not in keys:
                save.write(Line)
    save.close()
####


main()
