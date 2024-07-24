#!/usr/bin/env python3
import os
import sys
import argparse
import multiprocessing
import numpy as np
import pandas as pd
from functools import reduce
from operator import and_
from sklearn.manifold import TSNE
from sklearn import preprocessing
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import KNeighborsClassifier 
from sklearn.model_selection import train_test_split 
from sklearn.datasets import load_iris
from operator import itemgetter
import matplotlib.pyplot as plt
import pickle
parser = argparse.ArgumentParser(description='A Variety Recognition Model Based on Whole Genome SSR Digital Features')
parser.add_argument('-o','--output_path',required=True)
parser.add_argument('-s','--start_stage',required=True,type=int,choices=[1,2,3],default=1)
parser.add_argument('-e','--end_stage',required=True,type=int,choices=[1,2,3],default=3)
parser.add_argument('-miniap','--minia_proceedings',help='The process using while runing minia',required=False,type=int,default=8)
parser.add_argument('-miniac','--minia_cores',help='The core using while runing minia',required=False,type=int,default=10)
parser.add_argument('-spm','--spm_type',help='sequence_file_type,including single,pair,and mate-pair',required=False,choices=['single','pair','mate'],type=str,default='pair')
parser.add_argument('-misap','--misa_proceedings',help='The process using while runing misa',required=False,type=int,default=10)
parser.add_argument('-pp','--polymorphism_parameter',help='Parameters for calculating SSR polymorphism',required=True,type=float,default=0.8)
parser.add_argument('-index','--individual_variety_file',help='files that contains the relationship between individual(base_name) and variety',required=True)
args = vars(parser.parse_args())
output_path = os.path.abspath(args['output_path'])
start_stage = args['start_stage']
end_stage = args['end_stage']
minia_proceedings= args['minia_proceedings']
minia_cores=args['minia_cores']
misa_proceedings=args['misa_proceedings']
polymorphism_parameter=args['polymorphism_parameter']
spm_type=args['spm_type']
individual_variety_file=os.path.abspath(args['individual_variety_file'])
storage_name=os.path.basename(individual_variety_file)
print(storage_name)
print('minia_proceedings',minia_proceedings)
def get_individual_variety_dict(file_path):
    individual_variety_dict=dict()
    index_variety_dict=dict()
    all_individual=[]
    variety_list=[]
    with open(file_path) as file:
        index=0
        for line in file:
            line0=line.rstrip().split('\t')
            variety=line0[0]#get variety
            individual=line0[1]#get individual name
            index_variety_dict[index]=variety
            all_individual.append(individual)
            variety_list.append(variety)
            index+=1
            if variety not in individual_variety_dict:
                individual_variety_dict[variety]=[individual]
            else:
                individual_variety_dict[variety].append(individual)
    return individual_variety_dict,all_individual,variety_list,index_variety_dict
individual_variety_dict,all_individual,variety_list,index_variety_dict=get_individual_variety_dict(individual_variety_file)
varieties_number=int(len(set(variety_list)))
if start_stage > end_stage:
    print('error!start_stage must smaller than end stage')
    exit(1)
#the stage of operation
all_stage=[1,2,3]
opertional_stage=[]
for a in all_stage:
    if a >= start_stage and a <= end_stage:
        opertional_stage.append(a)
#build logfile,will overwrite earlier one
command = 'touch log.txt'
os.system(command)
command = 'echo log.txt > log.txt'
os.system(command)
#stage 1 minia
def mini(file_path,out):
    command = 'minia -in {} -kmer-size 31 -out-dir {} -nb-cores 10 -nb-glue-partitions 200'.format(file_path,out)
    print(command)
    os.system(command)
#stage 1-1 build inputfile of minia
def traverse_folder_build_inputfile(path):
    sequence_file_in_path = []
    sequence_file_prefix_list= []
    input_path_list = []
    for root, dirs, files in os.walk(path):
        for file in files:
            file0 = file.split('.',1)
            if len(file0) > 1:
                if file.split('.',1)[1] == 'fq':
                    sequence_file_in_path.append(file)
    sequence_file_in_path = list(set(sequence_file_in_path))
    for item in sequence_file_in_path:
        item0 = item.split('_R',1)[0]
        if item0 not in sequence_file_prefix_list:
            sequence_file_prefix_list.append(item0)
    for item0 in sequence_file_prefix_list:
        z = item0
        file1 = z + '_R1.fq'
        file2 = z + '_R2.fq'
        command1 = 'echo '+path+'/'+file1+'>'+path+'/'+item0
        command2 = 'echo '+path+'/'+file2+'>>'+path+'/'+item0
        os.system(command1)
        os.system(command2)
        input_path_list.append(path+'/'+item0)
    return sequence_file_prefix_list,input_path_list
if 1 in opertional_stage:
    if spm_type == 'pair':
        now_path=os.path.abspath('.')
        file_number=len(all_individual)
        log='totaly ' +str(file_number)+' sequence file'
        command = 'echo '+ log+'>> log.txt'
        os.system(command)
        pool = multiprocessing.Pool(minia_proceedings)
        for file in all_individual:
            intact_path = file
            log='minia is assemble '+ intact_path
            command = 'echo '+ log+'>> log.txt'
            os.system(command)
            pool.apply_async(mini, (intact_path,))
        pool.close()
        pool.join()
        for file_prefix in all_individual:
            command ='rm '+now_path+'/'+file_prefix+'.unitigs.fa'
            print(command)
            os.system(command)
            command ='rm '+file_prefix+'.h5'
            print(command)
            os.system(command)
            if os.path.abspath('.') != input_path:
                command ='mv '+now_path+'/'+file_prefix+'.contigs.fa ' +input_path
                print(command)
                os.system(command)
        log ='all unitigs file and h5 file have been deleted,and the contigs file is under input_path'
        print(log)
        command = 'echo '+log+'>> log.txt'
        os.system(command)
        log ='stage 1 has been ended'
        command = 'echo '+log+'>> log.txt'
        os.system(command)
#stage 2 misa
def misa(file):
    command ='misa.pl ' + file
    os.system(command)
    log ='misa.pl ' + file
    command = 'echo '+log+'>> log.txt'
    os.system(command)
def traverse_folder_contigs_file(path):
    contigs_file_path_list=[]
    for root, dirs, files in os.walk(path):
        for file in files:
            file0 = file.split('.',1)
            if len(file0) > 1:
                if file.split('.',1)[1] == 'contigs.fa':
                    print(path+'/'+file)
                    contigs_file_path_list.append(path+'/'+file)
    return contigs_file_path_list
if 2 in opertional_stage:
    pool = multiprocessing.Pool(40)
    for individual in all_individual:
        intact_path=individual
        pool.apply_async(misa, (intact_path, ))
    pool.close()
    pool.join()
    log ='misa has extracted all ssr,the misa_out file is in the input_path'
    command = 'echo '+log+'>> log.txt'
    os.system(command)
    log = 'stage 2 has been ended'
    command = 'echo '+log+'>> log.txt'
    os.system(command)
#stage 3 ssr_selection
def traverse_folder_misafile_path(path):
    lista=[]
    listb=[]
    for root, dirs, files in os.walk(path):
        for file in files:
            file0_split = file.split('.',1)
            if len(file0_split) > 1:
                if file.split('.',1)[1] == 'contigs.fa.misa':
                    lista.append(path+'/'+file)
                    listb.append(file)
    return lista,listb
def get_ssr_from_file(misafile_path):
    ssr_dict = dict()
    ssr_dict2 = dict()
    ssr_list=[]
    with open(misafile_path) as file:
        for line in file:
            if line.count ('(') == 1:
                r_part = line.split('(',1)[1]
                motifs = r_part.split(')',1)[0]
                motifs_hh = r_part.split(')',1)[1]
                motifs_n = motifs_hh.split('\t',1)[0]
                intact_ssr = '('+ motifs + ')'+motifs_n
                ssr_list.append(intact_ssr)
    ssr_number=0
    for ssr in ssr_list:
        ssr_number=ssr_number+1
        if ssr not in ssr_dict:
            ssr_dict[ssr]=1
        #else:
            #ssr_dict[ssr]+=1
    for ssr in ssr_dict:
        ssr_dict2[ssr]=(ssr_dict[ssr])/(ssr_number)
    return ssr_dict
def set_union(setA,setB):
    set_union = setA | setB
    return set_union
def standard_deviation(ListA):
    standard_deviation = np.std(ListA,ddof=1)
    return standard_deviation
def list_remove(listb,s):
    list_new=[]
    count=0
    for value in listb:
        if value==s:
            if count==0:
                count +=1
            else:
                list_new.append(value)
        else:
            list_new.append(value)
    return list_new
def ssr_selection(individual_variety_file):
    global storage_name
    individual_variety_dict,all_individual,variety_list,index_variety_dict=get_individual_variety_dict(individual_variety_file)
    directory_path=output_path+'/'+storage_name+'_'+str(polymorphism_parameter)
    print('directory_path',directory_path)
    try:
        os.mkdir(directory_path)
    except:
        pass
    log = 'now process the misa_out file'
    command = 'echo '+log+'>> log.txt'
    os.system(command)
    all_misa_file_path=[]
    for file in all_individual:
        file_misa=file+'.contigs.fa.misa'
        print(file_misa,'file_misa')
        all_misa_file_path.append(file_misa)
    variety_dict=dict()
    for variety in individual_variety_dict:
        for individual in individual_variety_dict[variety]:
            individual_path=individual+'.contigs.fa.misa'
            if variety not in variety_dict:
                variety_dict[variety]=[individual_path]
            else:
                variety_dict[variety].append(individual_path)
    nested_dict_include_ssr_info = dict()
    for misa_file in all_misa_file_path:
        nested_dict_include_ssr_info[misa_file]=get_ssr_from_file(misa_file)
    A =[set(nested_dict_include_ssr_info.get(individual)) for individual in all_misa_file_path]
    indivduals_union_ssr =list(reduce(set_union,A))
    all_ssr_data_base_on_union_ssr = dict()
    def get_occurrence_number_of_each_ssr_in_every_indivdual(SSR,individuals_path):
        for individual in individuals_path:
            if all_ssr_data_base_on_union_ssr.get(individual)==None:
                if nested_dict_include_ssr_info[individual].get(SSR)==None:
                    all_ssr_data_base_on_union_ssr[individual]=[0]
                else:
                    all_ssr_data_base_on_union_ssr[individual]=[nested_dict_include_ssr_info[individual].get(SSR)]
            else:
                if nested_dict_include_ssr_info[individual].get(SSR)==None:
                    all_ssr_data_base_on_union_ssr[individual]+=[0]
                else:
                    all_ssr_data_base_on_union_ssr[individual]+=[nested_dict_include_ssr_info[individual].get(SSR)]
    for ssr in indivduals_union_ssr:
        get_occurrence_number_of_each_ssr_in_every_indivdual(ssr,all_misa_file_path)
    all_ssr_data_base_on_union_ssr['SSR']= indivduals_union_ssr
    df = pd.DataFrame.from_dict(all_ssr_data_base_on_union_ssr)
    df.to_csv(output_path+'/all_data.csv',sep='\t',index=False)
    def calculate_the_polymorphism_of_ssr(ssr):
        count = 0
        for variety in variety_dict:
            M=[]
            for individual in all_misa_file_path:
                if nested_dict_include_ssr_info[individual].get(ssr)==None:
                    M.append(0)
                else:
                    M.append(nested_dict_include_ssr_info[individual].get(ssr))
            varietyA = []
            for individual in variety_dict.get(variety):
                if nested_dict_include_ssr_info[individual].get(ssr)==None:
                    varietyA.append(0)
                else:
                    varietyA.append(nested_dict_include_ssr_info[individual].get(ssr))
            variety_std = np.std(varietyA,ddof=1)
            for a in varietyA:
                M.remove(a)
            C_std = np.std(M,ddof=1)
            if variety_std < C_std:
                count +=1
        if (count/varieties_number)>polymorphism_parameter:
            return 1
    polymorphism_ssr=[]
    for ssr in indivduals_union_ssr:
        if calculate_the_polymorphism_of_ssr(ssr)==1:
            polymorphism_ssr.append(ssr)
    ########################################################################
    list_save_path=directory_path+'/'+'polymorphism_ssr_list.pkl'
    with open(list_save_path, 'wb') as file:
        pickle.dump(polymorphism_ssr,file)
    ########################################################################
    print('polymorphism_ssr',len(polymorphism_ssr))
    all_ssr_data_base_on_polymorphism_ssr=dict()
    all_ssr_data_base_on_polymorphism_ssr['SSR']=polymorphism_ssr
    def get_ssr_data_base_on_polymorphism_ssr(SSR,individuals_path):
        for individual in individuals_path:
            if all_ssr_data_base_on_polymorphism_ssr.get(individual)==None:
                if nested_dict_include_ssr_info[individual].get(SSR)==None:
                    all_ssr_data_base_on_polymorphism_ssr[individual]=[0]
                else:
                    all_ssr_data_base_on_polymorphism_ssr[individual]=[nested_dict_include_ssr_info[individual].get(SSR)]
            else:
                if nested_dict_include_ssr_info[individual].get(SSR)==None:
                    all_ssr_data_base_on_polymorphism_ssr[individual]+=[0]
                else:
                    all_ssr_data_base_on_polymorphism_ssr[individual]+=[nested_dict_include_ssr_info[individual].get(SSR)]
    pool = multiprocessing.Pool(10)
    for ssr in polymorphism_ssr:
        get_ssr_data_base_on_polymorphism_ssr(ssr,all_misa_file_path)
        pool.apply_async(get_ssr_data_base_on_polymorphism_ssr,(ssr,all_misa_file_path))
    pool.close()
    pool.join()
    df = pd.DataFrame.from_dict(all_ssr_data_base_on_polymorphism_ssr,orient='index')
    storage_name2='polymorphism_ssr'+str(polymorphism_parameter)+'.csv'
    storage_path=output_path+'/'+storage_name2
    df.to_csv(storage_path,sep='\t',index=True)
    def delete_first_line(filename):
        df = pd.read_csv(filename, skiprows=1)
        df.to_csv(filename, index=False)
    delete_first_line(storage_path)
    data=pd.read_csv(storage_path, sep = "\t")
    print(data)
    #############################Picture draw has down,model build strat#############################
    data=pd.read_csv(storage_path,sep='\t',index_col='SSR')
    X=data
    index=data.index
    #print('column',data.loc['/Rho/data/W-2-03.contigs.fa.misa'])
    y0=data.index
    y=variety_list
    X_train=X
    y_train=y
    estimator = KNeighborsClassifier(n_neighbors=7)
    estimator.fit(X,y)
    model_save_path =directory_path+'/'+'predict_model.pkl'
    print('path',model_save_path)
    with open(model_save_path, 'wb') as file:
        pickle.dump(estimator,file)
if 3 in opertional_stage:
    ssr_selection(individual_variety_file)
    '''
    individual_variety_dict2,all_individual2,variety_list2,index_variety_dict2=get_individual_variety_dict(all_individual_file)
    directory_path2=output_path+'/'+storage_name+'_'+str(polymorphism_parameter)
    for item in all_individual2:
        if item not in all_individual:
            full__path=item+'.contigs.fa.misa'
    command='python /Rho/pycat/Vimbwgs_model_predict.py -index '+individual_variety_file+' -i '+full__path+' -d '+directory_path2
    os.system(command)
    '''
