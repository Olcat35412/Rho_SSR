#!/usr/bin/env python3
import argparse
import pandas as pd
import csv
import os
from Bio import SeqIO
parser = argparse.ArgumentParser(description="extract genes habor ssr")
parser.add_argument('-txt','--txt_path',required=True)
parser.add_argument('-p','--pep_path',required=True)
parser.add_argument('-o','--output_path',required=True)
args = vars(parser.parse_args())
txt_path = os.path.abspath(args['txt_path'])
pep_path = os.path.abspath(args['pep_path'])
output_path =os.path.abspath(args['output_path'])
pep_name =os.path.basename(pep_path)
pep_name_no_ext=os.path.splitext(pep_name)[0]
txt_name=os.path.basename(args['txt_path'])
txt_name_no_ext=os.path.splitext(txt_name)[0]
with open (txt_path,'r') as file:
    column_list=[]
    for line in file:
        columns=line.rstrip('\n')
        column_list.append(columns)
print(column_list)
count=0
print(output_path)
with open(output_path,'w') as file:
    seq_records = SeqIO.parse(pep_path, "fasta")
    for seq_record in seq_records:
        #if seq_record.id.split(".")[0] in column_list:
        if seq_record.id in column_list:
            SeqIO.write(seq_record, file, "fasta")
            count += 1
print(count)
