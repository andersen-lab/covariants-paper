# script to get interesting cryptic evolution sequences
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import re
import json
from tqdm import tqdm

from outbreak_data import outbreak_data
from outbreak_tools import outbreak_tools

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

df = pd.read_csv('sd_cryptic.tsv',sep='\t')
import os
os.system('rm evoplots/*')
def getAAcoords(pos):
    return (pos - 21563)//3 +1

def sortFun_full(mut0):
    nuc0 = mut0.split('(S')[0].split('(O')[0]
    if '(' in nuc0:
        return int(nuc0.replace('(','').split(',')[0])+0.5
    else:
        return int(nuc0[1:len(nuc0)-1])

def sortFun(mut0):
    if 'DEL' in mut0:
        if '/' in mut0:
            return int(mut0.split('DEL')[1].split('/')[0])
        else:
            return int(mut0.split('DEL')[1]) 
    elif 'INS' in mut0:
        h0 = mut0.split('INS')[1]
        print(h0)
        return int(re.findall(r"(\d+)[A-Z]", h0)[0])
        # return int(mut0.split('INS')[1])
    else:
        # print(mut0)
        return int(mut0[3:(len(mut0)-1)])


def get_simplified_names(muts,dupMuts):
    #check if nonsynonymous

    # if type(muts) is not tuple:
    #     muts = tuple([muts,])
    mutsNew = []
    for mut0 in muts:
        if ')(' in mut0:
            mutsNew.append(mut0.split(')(')[1].split(')')[0])
        else:
            aa0 = mut0.split('(')[1].split(')')[0]
            part0 = aa0.split(':')[1]
            if part0[0]==part0[-1]:
                mutsNew.append(mut0.split('(')[0])
            else:
                mutsNew.append(aa0)
        if mutsNew[-1] in dupMuts:
            mutsNew[-1] = mut0
    return mutsNew

def get_duplicate_indices(lst):
    duplicates = {}
    for i, x in enumerate(lst):
        if x in duplicates:
            duplicates[x].append(i)
        else:
            duplicates[x] = [i]
    duplicates = [v for k,v in zip(duplicates.keys(),duplicates.values()) if len(v)>1]
    duplicates_flat = [d for d1 in duplicates for d in d1]
    return duplicates_flat

df['Covariants'] = df['Covariants'].apply(lambda x: tuple(set(x.replace("'","").replace('[','').replace(']','').split(', '))))
# df['Covariants_nuc'] = df['Covariants'].apply(lambda x: tuple([x0.split('(S')[0].split('(O')[0] for x0 in x if 'S:' in x0]))

#remove frame shifting deletions 
# df['Covariants_AA'] = df['Covariants'].apply(lambda x: tuple(['S' + x0.split('(S')[1].split(')')[0] for x0 in x if 'S:' in x0]))
df['collection_date'] = pd.to_datetime(df['collection_date'])
df = df[df['collection_date']<'2024-11-01']
# df['Coverage_start'] = df['Coverage_start'].apply(lambda x:getAAcoords(x))
# df['Coverage_end'] = df['Coverage_end'].apply(lambda x:getAAcoords(x))

df['num_muts'] = df['Covariants'].apply(lambda x:len(x))
df = df[df['num_muts']>=2]
df = df.sort_values(by='num_muts',ascending=False).reset_index(drop=True)

## re-sort to ensure we don't miss matches. 
# df['Covariants_AA'] = df['Covariants_AA'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun)))
df['Covariants'] = df['Covariants'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun_full)))

# df = df[df['Count']>=25]
df1 = df.groupby('Covariants').agg({'Count':tuple,'location':tuple,'collection_date':tuple,
                                       'Coverage_start':tuple,'Coverage_end':tuple,
                                       'num_clinical_detections':'mean','num_muts':'min'}).sort_values(by='num_muts',ascending=False)
df1['Total_Obs']= df1['Count'].apply(lambda x: len(x))


# pull only cryptics (here, defined as 1 or less clinical detections)

df_cryptics = df1[df1['num_clinical_detections']<=1]
# df_cryptics = df_cryptics[df_cryptics['num_muts']>=1]
df_cryptics['Covariants'] = df_cryptics.index
df_cryptics = df_cryptics[df_cryptics['Total_Obs']>=2]
# remove cryptics that are subsets of other cryptics. 
j=0
while j<df_cryptics.shape[0]:
    test_clust = set(df_cryptics.index[j])
    df_subset = df_cryptics[df_cryptics['Covariants'].apply(lambda x: set(x).issubset(test_clust) & (set(x)!=test_clust))]
    df_cryptics = df_cryptics.drop(index=df_subset.index)
    j=j+1

## now generate subtrees leading up to these cryptics. 

# df1 = df1[df1['num_clinical_detections'] >=1]
# count_dict = {}
for i0 in tqdm(range(0,df_cryptics.shape[0])):
    test_clust = set(df_cryptics.index[i0])
    #for testing, enforce complete of all mutations in the most mutated. 
    df0 = df.copy()
    positions = np.array([sortFun_full(t) for t in test_clust])
    df0 = df0[df0['Coverage_start']<=positions.min()]
    df0 = df0[df0['Coverage_end']>=positions.max()]
    # now group
    # df0 = df0[df0['Count']>=50]

    df0 = df0.groupby('Covariants').agg({'Count':tuple,'location':tuple,'collection_date':tuple,
                                            'Coverage_start':tuple,'Coverage_end':tuple,
                                            'num_clinical_detections':'mean',
                                            'num_muts':'min'})
                                            
    df0['Total_Obs']= df0['Count'].apply(lambda x: len(x))
    df0 = df0[df0['Total_Obs']>=2]
    df0['Covariants'] = df0.index
    # df0 = df0.reset_index(drop=False)
    df0['num_muts'] = df0['Covariants'].apply(lambda x:len(x))
    df0 = df0.sort_values(by='num_muts',ascending=False)

    df_superset = df0[df0['Covariants'].apply(lambda x: set(x).issubset(test_clust))]
    # print(i0,df_superset)
    # df_superset = df[df['Covariants'].apply(lambda x: test_clust.issubset(set(x)))]
    # df_superset = pd.concat((df_subset,df_superset),axis=0).drop_duplicates().sort_values(by='num_muts')
    if df_superset.shape[0]<=1:
        print(f'No detected stepwise evolution for {test_clust}')
        continue

    df_superset.loc[:,'Covariants'] = df_superset['Covariants'].apply(lambda x: tuple(sorted(list(x),key=sortFun_full)))

    hh = list(df_superset.index)
    hh_flat = list(set([h0 for hh0 in hh for h0 in hh0]))
    hh_simplified = get_simplified_names(hh_flat,[])
    if len(set(hh_simplified))< len(hh_flat):
        dupInds = get_duplicate_indices(hh_simplified)
        dupMuts = set([hh_simplified[i0] for i0 in dupInds])
    else:
        dupMuts = set([])

    parent_list = []
    new_muts_list = []
    for j,c in enumerate(df_superset['Covariants']):
        # get subsets
        subsets = df_superset[df_superset['Covariants'].apply(lambda x: set(x).issubset(c) & (x!=c))]
        if subsets.shape[0]>0:
            subsets = pd.concat((subsets,pd.Series(subsets['Covariants'].apply(lambda x: len(set(c)-set(x))),name='difference')),axis=1)
            minDiff = subsets['difference'].min()
            subsets = subsets[subsets['difference']==minDiff]
            
            if subsets.shape[0] ==1:
                parent_list.append([subsets['Covariants'].iloc[0]])
                new_muts_list.append([tuple(set(c)-set(subsets['Covariants'].iloc[0]))])
            else:
                parent_list.append(list(subsets['Covariants']))
                new_muts_list.append([tuple(set(c)- set(t0)) for t0 in subsets['Covariants']])
            # print(c,subsets)
        else:
            parent_list.append(None)
            new_muts_list.append(c)

    from networkx.drawing.nx_agraph import graphviz_layout
    # for each, define parent, unless smallest cluster. 
    edge_labels = {}
    G = nx.DiGraph()
    # assemble graph
    # make edges to first level
    for j,c in enumerate(df_superset['Covariants']):
        if parent_list[j] is not None:
            for l in range(0,len(parent_list[j])):
                G.add_edge(parent_list[j][l],c,weight=1)
                edge_labels[(parent_list[j][l],c)] = ','.join(get_simplified_names(new_muts_list[j][l],dupMuts))
    

    try:
        if len(nx.dag_longest_path(G))<=2:
            continue
    except:
        pass
    fig,ax = plt.subplots()#figsize=(5,1.5))
    pos=graphviz_layout(G, prog='dot',args="-Grankdir='LR' -Goverlap=false",root=0)
    nx.draw_networkx_nodes(G, pos, node_color="cornflowerblue",alpha=0.9,margins=0.1,node_size=80)
    nx.draw_networkx_edges(
        G,
        pos,
        width=1,
        alpha=0.8,
        arrows=True,
        edge_color="grey",
    )
    labels = {c:o for c,o in zip(df_superset['Covariants'],list(df_superset['Total_Obs']))}
    if len(labels)> len(pos):
        labels = {key:labels[key] for key in pos.keys()}
    nx.draw_networkx_labels(G, pos = pos, labels = labels, font_color="white", font_size=6)#,verticalalignment='bottom')
    nx.draw_networkx_edge_labels(G, pos = pos, edge_labels=edge_labels, font_color='black',font_size=5)
    minX = 500
    maxX = 900

    ymin, ymax = ax.get_ylim()
    y_length = ymax - ymin
    
    xmin, xmax = ax.get_xlim()
    x_length = xmax - xmin

    start = df_superset['Covariants'].iloc[-1]
    minXpos = np.min([x[0] for x in pos.values()])

    for pk in pos.keys():
        if G.in_degree(pk)==0:
            # ax.text(pos[pk][0],pos[pk][1]+y_length*0.04,'\n'.join(pk),fontsize=5,verticalalignment='bottom',horizontalalignment='center')
            # else:
            ax.text(pos[pk][0]-x_length*0.06,pos[pk][1],'\n'.join(get_simplified_names(pk,dupMuts)),fontsize=5,verticalalignment='center',horizontalalignment='right')
        if G.out_degree(pk)==0:
            # ax.text(pos[pk][0],pos[pk][1]+y_length*0.04,'\n'.join(pk),fontsize=5,verticalalignment='bottom',horizontalalignment='center')
            # else:
            ax.text(pos[pk][0]+x_length*0.06,pos[pk][1],'\n'.join(get_simplified_names(pk,dupMuts)),fontsize=5,verticalalignment='center',horizontalalignment='left')
    # add number of clinical detections below ww detection count
    clin_detects = {c:counts for c,counts in zip(df_superset['Covariants'],df_superset['num_clinical_detections'])}
    for pk in pos.keys():
        ax.text(pos[pk][0],pos[pk][1]-y_length*0.03,int(df0.loc[[pk],'num_clinical_detections']),
                fontsize=5,verticalalignment='center',horizontalalignment='center',color='red')
    # add the sites the variant was detected at. 
    for pk in pos.keys():
        locs = set(list(df0.loc[[pk],'location'][0]))
        for j,l0 in enumerate(['South Bay','Point Loma','Encina']):
            if l0 in locs:
                ax.scatter(pos[pk][0]+x_length*0.02,pos[pk][1]+ y_length*0.02*(1-j),marker='s',s=5,color='darkblue',linewidth=0.5)
            else:
                ax.scatter(pos[pk][0]+x_length*0.02,pos[pk][1]+ y_length*0.02*(1-j),marker='s',s=5,edgecolors='darkblue',facecolors='none',linewidth=0.5)
    
    for pk in pos.keys():
        tp = pd.Series(df0.loc[[pk],'collection_date'][0],name='times')
        if len(tp)==1:
                ax.text(pos[pk][0],pos[pk][1]+y_length*0.05,
                f"{tp.iloc[0].strftime('%Y-%m-%d')}",
                fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')
        else:
            ax.text(pos[pk][0],pos[pk][1]+y_length*0.05,
                    f"{tp.min().strftime('%Y-%m-%d')}\n-{tp.max().strftime('%Y-%m-%d')}",
                    fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')

    plt.gca().set_frame_on(False)

    plt.savefig(f'evoplots/ww_evo_seq{i0}.pdf',transparent=True)
    plt.close('all')

