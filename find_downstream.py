# # script to get interesting cryptic evolution sequences
# import pandas as pd
# import numpy as np
# import networkx as nx
# import matplotlib.pyplot as plt
# import matplotlib
# import re
# import json

# from outbreak_data import outbreak_data
# from outbreak_tools import outbreak_tools

# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42


# df = pd.read_csv('all_covariants_cryptic_new.csv')

# def getAAcoords(pos):
#     return (pos - 21563)//3 +1

# def sortFun(mut0):
#     if 'DEL' in mut0:
#         if '/' in mut0:
#             return int(mut0.split('DEL')[1].split('/')[0])
#         else:
#             return int(mut0.split('DEL')[1]) 
#     elif 'INS' in mut0:
#         h0 = mut0.split('INS')[1]
#         print(h0)
#         return int(re.findall(r"(\d+)[A-Z]", h0)[0])
#         # return int(mut0.split('INS')[1])
#     else:
#         return int(mut0[3:(len(mut0)-1)])
    
# def correct_del_name(mut0):
#     if 'DEL' in mut0:
#         if '/' in mut0:
#             return mut0
#         else:
#             return mut0 + '/' + mut0.split('DEL')[1]
#     else:
#         return mut0
    
# # df = df[df['Coverage_start']<=22810]
# # df = df[df['Coverage_end']>=22950]
# def find_frameshift(mut0):
#     if ')' in mut0:
#         mut0 = mut0[0:(len(mut0)-1)].split(',')[1]
#         if mut0[0]=="'":
#             mut0 = mut0.replace("'","")
#             if len(mut0)%3 ==0:
#                 return True
#             else:
#                 return False
#         else:
#             if int(mut0)%3 ==0:
#                 return True
#             else:
#                 return False              
#     else:
#         return True

# count_dict = {}
# with open('count_data.json', 'r') as file:
#     count_dict = json.load(file)
# # count_dict = {}
# df['Covariants'] = df['Covariants'].apply(lambda x: tuple(set(x.split(' '))))
# df['Covariants_nuc'] = df['Covariants'].apply(lambda x: tuple([x0.split('(S')[0].split('(O')[0] for x0 in x]))
# #remove frame shifting deletions 
# df['Covariants_AA'] = df['Covariants'].apply(lambda x: tuple(['S' + x0.split('(S')[1].split(')')[0] if 'S:' in x0 else '' for x0 in x]))
# df['keep_inds'] = df['Covariants_nuc'].apply(lambda x: [i for i,x0 in enumerate(x) if find_frameshift(x0)])

# df['Covariants_nuc'] = df[['keep_inds','Covariants_nuc']].apply(lambda x: [x['Covariants_nuc'][i0] for i0 in x['keep_inds'] if len(x['Covariants_nuc'][i0])>0],axis=1)
# df['Covariants_AA'] = df[['keep_inds','Covariants_AA']].apply(lambda x: tuple([x['Covariants_AA'][i0] for i0 in x['keep_inds'] if len(x['Covariants_AA'][i0])>0]),axis=1)
# df['Covariants'] = df[['keep_inds','Covariants']].apply(lambda x: [x['Covariants'][i0] for i0 in x['keep_inds']],axis=1)

# df['Covariants'] = df['Covariants_AA']
# # df['Covariants']  = df['Covariants'].apply(lambda x: tuple(set(x.replace("[",'').replace("]",'').replace("'","").split(', '))))
# df['Coverage_start'] = df['Coverage_start'].apply(lambda x:getAAcoords(x))
# df['Coverage_end'] = df['Coverage_end'].apply(lambda x:getAAcoords(x))


# # df = df.groupby('Covariants').agg({'Count':tuple,'location':tuple,'collection_date':tuple,'Coverage_start':tuple,'Coverage_end':tuple,'num_clinical_detections':'mean'})
# # df0['Total_Obs']= df0['Count'].apply(lambda x: len(x))
# # df0 = df0.reset_index(drop=False)

# df['num_muts'] = df['Covariants'].apply(lambda x:len(x))
# # df = df[df['num_muts']>=2]
# # df = df[df['Count']>10]
# df = df.sort_values(by='num_muts',ascending=False).reset_index(drop=True)
# df['Covariants'] = df['Covariants'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun)))
# df = df.groupby('Covariants').agg({'Count':tuple,'location':tuple,'collection_date':tuple,'Coverage_start':tuple,'Coverage_end':tuple,'num_clinical_detections':'mean'}).reset_index()
# df['Total_Obs']= df['Count'].apply(lambda x: len(x))

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
from os import system
system('rm descent_plots/*.pdf')
df = pd.read_csv('sd_cryptic.tsv',sep='\t')
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

df = df.drop_duplicates()
df['Covariants'] = df['Covariants'].apply(lambda x: tuple(set(x.replace("'","").replace('[','').replace(']','').split(', '))))
# df['Covariants_nuc'] = df['Covariants'].apply(lambda x: tuple([x0.split('(S')[0].split('(O')[0] for x0 in x if 'S:' in x0]))

#remove frame shifting deletions 
# df['Covariants_AA'] = df['Covariants'].apply(lambda x: tuple(['S' + x0.split('(S')[1].split(')')[0] for x0 in x if 'S:' in x0]))
df['collection_date'] = pd.to_datetime(df['collection_date'])
# df['Coverage_start'] = df['Coverage_start'].apply(lambda x:getAAcoords(x))
# df['Coverage_end'] = df['Coverage_end'].apply(lambda x:getAAcoords(x))
df['Covariants'] = df['Covariants'].apply(lambda x: tuple(sorted(list(set(x)),key=sortFun_full)))
df['Covariants'] = df['Covariants'].apply(lambda x:tuple(get_simplified_names(x,[])))

df['num_muts'] = df['Covariants'].apply(lambda x:len(x))
# df = df[df['num_muts']>=2]
df = df[df['Count']>=50]
df = df.sort_values(by='num_muts',ascending=False).reset_index(drop=True)
# hh = list(df['Covariants'])
# from collections import Counter
# print(Counter([h0 for hh0 in hh for h0 in hh0]))

df = df.groupby('Covariants').agg({'Count':tuple,'location':tuple,'collection_date':tuple,'Coverage_start':tuple,'Coverage_end':tuple,'num_clinical_detections':'mean'})
df['Covariants'] = df.index
df['Total_Obs']= df['Count'].apply(lambda x: len(x))
df = df[df['Total_Obs']>=2]
## for a known mutation cluster, see what evolves from it. 
clusters = [("S:G142D","S:DEL144","S:F157S","S:R158G"),("S:H655Y", "S:N679K", "S:P681R"),('S:K356T', 'S:S371F', 'S:S373P', 'S:S375F',"S:D405N","S:R408S"),('S:N969K','S:Q954H')]#("S:S371F","S:S373P","S:S375F")#("S:K356T","S:S371F","S:S373P","S:S375F","S:T376A","S:L390F","S:R403K")#('S:G142D', 'S:DEL144', 'S:F157S', 'S:R158G')
for cluster0 in clusters:
    #,"S:DEL463") #("S:K356T","S:S371F","S:S373P","S:S375F","S:T376A","S:L390F","S:R403K")#("S:K417N","S:N440K","S:V445P")#("S:H655Y","S:N679K","S:P681H","S:A653T")
                # ("S:G142D","S:DEL144","S:F157S","S:R158G")#,"S:T167I") #"S:C136F"
    #("S:H655Y","S:N679K","S:P681H","S:A653T") #("S:S477N","S:T478K","S:DEL483","S:E484K","S:F486P","S:Q498R","S:N501Y") #("S:H655Y","S:T678A","S:N679K","S:P681R") #("S:H655Y","S:N679K","S:P681H","S:A653T") ##
    df_superset = df[df['Covariants'].apply(lambda x: set(cluster0).issubset(set(x)) or set(cluster0)==set(x))]

    df_superset['num_descendants'] = [df_superset[df_superset['Covariants'].apply(lambda x: set(c0).issubset(set(x)))].shape[0] for c0 in df_superset.index] 
    df_superset = df_superset.sort_values(by='Total_Obs',ascending=False)


    if df_superset.shape[0]>10:
        df_superset = df_superset.iloc[0:10]

    #force the seed cluster to be in the matrix
    if df_superset[df['Covariants'].apply(lambda x: set(cluster0)==set(x))].shape[0]==0:
        df_superset = df_superset.reset_index(drop=True)
        newRow = pd.DataFrame({c0:None for c0 in df.columns},index=[cluster0])
        newRow['Covariants'] =[cluster0]
        df_superset = pd.concat((df_superset,newRow),axis=0,ignore_index=True).set_index('Covariants')
        df_superset['Covariants'] = df_superset.index

    # df_superset = df_superset[df_superset['num_clinical_detections']<=100 | (df_superset['Covariants']==cluster0)]
    parent_list = []
    new_muts_list = []
    for j,c in enumerate(df_superset['Covariants']):
        # get subsets
        subsets = df_superset[df_superset['Covariants'].apply(lambda x: set(x).issubset(c) & (x!=c))]
        if subsets.shape[0]>0:
            subsets['difference'] = subsets['Covariants'].apply(lambda x: len(set(c)-set(x)))
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
                edge_labels[(parent_list[j][l],c)] = ','.join(new_muts_list[j][l])


    l0 = len(nx.dag_longest_path(G))-1
    fig,ax = plt.subplots(figsize=(3*l0,5))
    # pos = nx.spring_layout(G, seed=1)
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

    labels.pop(cluster0)
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
            ax.text(pos[pk][0]-x_length*0.03,pos[pk][1],'\n'.join(pk),fontsize=5,verticalalignment='center',horizontalalignment='right')
    #     if G.out_degree(pk)==0:
    #         ax.text(pos[pk][0]+x_length*0.06,pos[pk][1],'\n'.join(pk),fontsize=5,verticalalignment='center',horizontalalignment='left')

    for pk in pos.keys():
        if df_superset.loc[[pk],'location'][0] is not None:
            locs = set(list(df_superset.loc[[pk],'location'][0]))
        else:
            continue
        for j,l0 in enumerate(['South Bay','Point Loma','Encina']):
            if l0 in locs:
                ax.scatter(pos[pk][0]+x_length*0.035,pos[pk][1]+ y_length*0.015*(1-j),marker='s',s=5,color='darkblue',linewidth=0.5)
            else:
                ax.scatter(pos[pk][0]+x_length*0.035,pos[pk][1]+ y_length*0.015*(1-j),marker='s',s=5,edgecolors='darkblue',facecolors='none',linewidth=0.5)
    # add number of clinical detections below ww detection count
    clin_detects = {c:counts for c,counts in zip(df_superset['Covariants'],df_superset['num_clinical_detections'])}
    for pk in pos.keys():
        ax.text(pos[pk][0],pos[pk][1]-y_length*0.03,str(int(df_superset.loc[[pk],'num_clinical_detections'])) if pk!=cluster0 else "",
                fontsize=5,verticalalignment='center',horizontalalignment='center',color='red')
        
    for pk in pos.keys():
        if df_superset.loc[[pk],'location'][0] is not None:
            tp = pd.Series(df_superset.loc[[pk],'collection_date'][0],name='times')
        else:
            continue
        if len(tp)==1:
                ax.text(pos[pk][0],pos[pk][1]+y_length*0.03,
                f"{tp.iloc[0].strftime('%Y/%m/%d')}",
                fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')
        else:
            ax.text(pos[pk][0],pos[pk][1]+y_length*0.03,
                    f"{tp.min().strftime('%Y/%m/%d')}\n-{tp.max().strftime('%Y/%m/%d')}",
                    fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')
            
    plt.gca().set_frame_on(False)
    fn0 = ','.join(cluster0)
    fig.tight_layout()
    plt.savefig(f'descent_plots/ww_evo_seq{fn0}.pdf',transparent=True)
    plt.close('all')