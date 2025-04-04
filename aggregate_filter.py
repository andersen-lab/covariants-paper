# script to get interesting cryptic evolution sequences
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import re
import json
from tqdm import tqdm
import os

from outbreak_data import outbreak_data
from outbreak_tools import outbreak_tools

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def getAAcoords(pos):
    # just for spike
    return (pos - 21563)//3 +1

def parse_cov_read(dir_,fn0):
    df = pd.read_csv(os.path.join(dir_, fn0), sep='\t')
    df['Collection_date'] = fn0.split('.')[0].split('_')[0]
    df['Location'] = fn0.split('.')[0].split('_')[1]
    return df

def find_frameshift(mut0):
    if ')' in mut0:
        mut0 = mut0[0:(len(mut0)-1)].split(',')[1]
        if mut0[0]=="'":
            mut0 = mut0.replace("'","")
            if len(mut0)%3 ==0:
                return True
            else:
                return False
        else:
            if int(mut0)%3 ==0:
                return True
            else:
                return False              
    else:
        return True


def sortFun(mut0):
    if 'DEL' in mut0:
        if '/' in mut0:
            return int(mut0.split('DEL')[1].split('/')[0])
        else:
            return int(mut0.split('DEL')[1]) 
    elif 'INS' in mut0:
        h0 = mut0.split('INS')[1]
        # print(h0)
        return int(re.findall(r"(\d+)[A-Z]", h0)[0])
        # return int(mut0.split('INS')[1])
    else:
        return int(mut0[3:(len(mut0)-1)])

def sortFun_nuc(x):
    # sort based on nuc position, ignoring nuc identities
    if ')' in x:
        return int(x[1:(len(x)-1)].split(',')[0])           
    else:
        return int(x[1:(len(x)-1)])

def correct_del_name(mut0):
    if 'DEL' in mut0:
        if '/' in mut0:
            return mut0
        else:
            return mut0 + '/' + mut0.split('DEL')[1]
    else:
        return mut0

print('Loading')
dir0 = 'SD-Wastewater-Cryptic-Variants/covariants'
df = pd.concat([parse_cov_read(dir0,fn) for fn in os.listdir(dir0) if fn.startswith('2')], ignore_index=True)
print('Cleaning up')
# count_dict = {}
df['Covariants'] = df['Covariants'].apply(lambda x: tuple(set(x.split(' '))))
df['Covariants_nuc'] = df['Covariants'].apply(lambda x: tuple([x0.split('(S')[0].split('(O')[0] for x0 in x]))
#remove frame shifting deletions 
df['Covariants_AA'] = df['Covariants'].apply(lambda x: tuple(['S' + x0.split('(S')[1].split(')')[0] if 'S:' in x0 else '' for x0 in x]))
df['keep_inds'] = df['Covariants_nuc'].apply(lambda x: [i for i,x0 in enumerate(x) if find_frameshift(x0)])


df['Covariants_nuc'] = df[['keep_inds','Covariants_nuc']].apply(lambda x: [x['Covariants_nuc'][i0] for i0 in x['keep_inds'] if len(x['Covariants_nuc'][i0])>0],axis=1)
df['Covariants_AA'] = df[['keep_inds','Covariants_AA']].apply(lambda x: tuple([x['Covariants_AA'][i0] for i0 in x['keep_inds'] if len(x['Covariants_AA'][i0])>0]),axis=1)
df['Covariants'] = df[['keep_inds','Covariants']].apply(lambda x: [x['Covariants'][i0] for i0 in x['keep_inds']],axis=1)

df = df[df['Covariants'].apply(lambda x: len(x)>0)]

#sort and merge back
df['Covariants_nuc'] = df['Covariants_nuc'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun_nuc)))
df['Covariants_AA'] = df['Covariants_AA'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun)))
# df['Covariants'] = df[['Covariants_nuc','Covariants_AA']].apply(lambda x:t[f'{xn}({xa})' for xn,xa in zip(x['Covariants_nuc'],x['Covariants_AA'])])
# df['Covariants']  = df['Covariants'].apply(lambda x: tuple(set(x.replace("[",'').replace("]",'').replace("'","").split(', '))))
df['Coverage_start'] = df['Coverage_start'].apply(lambda x:getAAcoords(x))
df['Coverage_end'] = df['Coverage_end'].apply(lambda x:getAAcoords(x))

## remove empty rows (entirely composed of frameshifted indels, and thus probably not real)
df['num_muts'] = df['Covariants'].apply(lambda x:len(x))
df = df[df['num_muts']>=2]
df = df.sort_values(by='num_muts',ascending=False).reset_index(drop=True)
df = df[df['Count']>=25]

print('Grouping')
df1 = df.groupby(['Covariants_nuc','Covariants_AA']
                 ).agg({'Count':tuple,'Location':tuple,
                        'Collection_date':tuple,'Coverage_start':tuple,
                        'Coverage_end':tuple,
                        'num_muts':'mean'}).sort_values(by='num_muts'
                                                        ,ascending=False).reset_index(drop=False)

# pull definitions from
# https://raw.githubusercontent.com/andersen-lab/Freyja/refs/heads/main/freyja/data/lineage_mutations.json

with open('lineage_mutations.json', 'r') as file:
    lindefs = json.load(file)

# cluster0 = ['C28687A','G28881A','G28882A','G28883C']
# lind = lindefs['AJ.1']
screen = False
if screen:
    def check_if_existing(cluster,lindef0):
        lindef = [l for l in lindef0 if ',' not in l]
        if cluster[0] in lindef:
            startInd = lindef.index(cluster[0])
            if startInd + len(cluster) < len(lind):
                lin0 = lindef[startInd:startInd+len(cluster)]
                if lin0 == list(cluster):
                    return True
        return False

    print('running screen')
    # check if cluster is a likely part of a known lineage
    check_list = [False]*df1.shape[0]
    for j,clust in tqdm(enumerate(df1['Covariants_nuc'])):
        for lind in lindefs.values():
            if check_if_existing(clust,lind):
                check_list[j]=True
                break
    df1['in_barcodes'] = check_list

hit_count = []
for mut0 in tqdm(df1['Covariants_AA']):
    pk0 = [correct_del_name(p0) for p0 in list(mut0)]
    hits = outbreak_data.mutation_prevalences(','.join(pk0))
    if hits.shape[0]>0:
        total_hits = hits['mutation_count'].sum()
    else:
        total_hits = 0
    hit_count.append(total_hits)

df1['clinical_detections'] = hit_count
df1.to_csv('cluster_counts.tsv',sep='\t')