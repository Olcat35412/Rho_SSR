#!/usr/bin/env python3
import argparse
import os
parser = argparse.ArgumentParser(description='cut_off_ssr_info')
parser.add_argument('-i','--input_path',required=True)
parser.add_argument('-o','--output_path',required=True)
parser.add_argument('-ssr','--ssr_condition',required=True)
parser.add_argument('-cds','--cds_condition',required=True)
args = vars(parser.parse_args())
input_path = os.path.abspath(args['input_path'])
output_path = os.path.abspath(args['output_path'])
med_path=output_path+".csv"
ssr_condition=int(args['ssr_condition'])
cds_condition=int(args['cds_condition'])
def filter_table(input_file,output_file,ssr_condition,cds_condition):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = next(f_in)
        f_out.write(header)
        for line in f_in:
            columns = line.strip().split('\t')
            ssr_value = int(columns[2])
            cds_value = int(columns[5])
            if ssr_value >= ssr_condition and cds_value >= cds_condition:
                f_out.write(line)
filter_table(input_path,med_path,ssr_condition,cds_condition)
command="awk -F'\t' 'NR>1 {print $1}' "+med_path+">"+output_path
os.system(command)
