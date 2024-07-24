#!/usr/bin/env python3
import argparse
import pandas as pd
import csv
import os
from Bio import SeqIO
parser = argparse.ArgumentParser(description="annotate CDS with SSR and pep translated form SSR")
parser.add_argument('-cds','--cds_path',required=True)
parser.add_argument('-ssr','--ssr_path',required=True)
parser.add_argument('-o','--output_path',required=True)
args = vars(parser.parse_args())
cds_path = os.path.abspath(args['cds_path'])
ssr_path = os.path.abspath(args['ssr_path'])
output_path = os.path.abspath(args['output_path'])
ssr_data = pd.read_csv(ssr_path, sep='\t')
ssr_info = list(ssr_data.iloc[:, [4,5,6,0]].itertuples(index=False, name=None))
def extract_sequence_from_fasta(fasta_file, seq_id, start, end):
    sequence = ""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == seq_id:
            sequence = record.seq[start-1:end]
            print(sequence)
            break
    protein_sequence = sequence.translate()
    return protein_sequence
def get_sequence_length(fasta_file, sequence_id):
    records = SeqIO.parse(fasta_file, "fasta")
    for record in records:
         if record.id == sequence_id:
             return len(record.seq)
def find_nearest_multiple_of_three_interval(a, b):
    nearest_a = a + ((3 - (a % 3) + 1) % 3)
    nearest_b = b - (b % 3) if b % 3 != 0 else b
    c = nearest_a if nearest_a >= a else nearest_b
    d = nearest_b if nearest_b <= b else nearest_a
    return c,d
class ssr:
    def __init__(self,cds_id,start_site,end_site,ssr_length,ssr_sequence):
        self.cds_id=cds_id
        self.start_site=start_site
        self.end_site=end_site
        self.ssr_length=ssr_length
        self.ssr_sequence=ssr_sequence
    @property
    def pep(self):
        new_start,new_end=find_nearest_multiple_of_three_interval(int(self.start_site),int(self.end_site))
        print(new_start,new_end,self.cds_id)
        pep=extract_sequence_from_fasta(cds_path,self.cds_id,new_start,new_end)
        return pep
    @property
    def cds_length(self):
        cds_length=get_sequence_length(cds_path,self.cds_id)
        return cds_length
dict_CDS_SSR_ANNO=dict()
count=0
with open(ssr_path) as file:
    next(file)
    for line in file:
        count+=1
        lines = line.rstrip('\n').split('\t')
        ssr_instance=ssr(lines[0].split('_')[0],lines[5],lines[6],lines[4],lines[3])
        if ssr_instance.cds_id not in dict_CDS_SSR_ANNO:
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]={}
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_site']=[]
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_length']=0
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_sequence']=[]
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_pep']=[]
            dict_CDS_SSR_ANNO[ssr_instance.cds_id]['cds_length']=ssr_instance.cds_length
        dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_site'].append((int(ssr_instance.start_site),int(ssr_instance.end_site)))
        dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_length']+=int(ssr_instance.ssr_length)
        dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_pep'].append(ssr_instance.pep)
        dict_CDS_SSR_ANNO[ssr_instance.cds_id]['ssr_sequence'].append(ssr_instance.ssr_sequence)
        #if count==100:
            #break
df=pd.DataFrame(dict_CDS_SSR_ANNO).T
df.to_csv(output_path,index_label='key',sep='\t')
print(output_path)
