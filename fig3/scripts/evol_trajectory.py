import pandas as pd
import numpy as np
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

import matplotlib.pyplot as plt
import matplotlib

""" Given a cryptic mutation cluster, plot the stepwise evolution of the mutations that descent from it."""

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def get_duplicate_indices(indices):
    duplicates = {}
    for i, x in enumerate(indices):
        if x in duplicates:
            duplicates[x].append(i)
        else:
            duplicates[x] = [i]
    duplicates = [v for k,v in zip(duplicates.keys(),duplicates.values()) if len(v)>1]
    duplicates_flat = [d for d1 in duplicates for d in d1]
    return duplicates_flat

def parse_query_list(query_list):
    query_list = query_list.strip('[]').split(', ')
    query_list = [x.strip().strip("'").strip('"') for x in query_list]
    return query_list


def get_aa_site(mut):
    if 'DEL' in mut:
        return int(mut.split('DEL')[1].split('/')[0])
    else: 
        return int(mut.split(':')[1][1:-1]
)

def main():
    df = pd.read_csv('../../data/cryptic_variants_metadata.tsv',sep='\t')

    df['Covariants'] = df['query'].apply(parse_query_list)
    df['Covariants'] = df['Covariants'].apply(lambda x: tuple(sorted(x,key=get_aa_site)))
    df['collection_date'] = pd.to_datetime(df['collection_date'])

    df['num_muts'] = df['Covariants'].apply(lambda x:len(x))
    df = df[df['num_muts']>=2]
    df = df.sort_values(by='num_muts',ascending=False).reset_index(drop=True)

    df = df[df['count']>=10]

    df.drop_duplicates(subset=['Covariants', 'collection_date', 'location'],keep='first',inplace=True)

    # convert coverage start and end to amino acid positions
    df['coverage_start'] = df['coverage_start'].apply(lambda x: (x - 21563) // 3 )
    df['coverage_end'] = df['coverage_end'].apply(lambda x: (x - 21563) // 3 )


    # aggregate metadata for each unique cluster
    df = df.groupby('Covariants').agg({'count':tuple,'location':tuple,'collection_date':tuple,'coverage_start':tuple,'coverage_end':tuple,'num_clinical_detections':'mean'})
    
    df['Covariants'] = df.index
    df['Total_Obs']= df['count'].apply(lambda x: len(x))
    df = df[df['Total_Obs']>=2]

    df.to_csv('df_filtered.tsv',sep='\t',index=False)

    ## for a known mutation cluster, see what evolves from it. 
    clusters = [
                #('S:G142D', 'S:DEL144/144', 'S:H146Q'),
                ("S:D614G", "S:P621S", "S:H655Y")
                #("S:H655Y", "S:N679K")
                # ("S:W633R"),
                # ("S:N679K", "S:S698P"),
                # ("S:M697V", "S:N679K", "S:P681H"),
                # ("S:H655Y","S:T678A","S:N679K","S:P681R"),
                # ('S:K356T', 'S:S371F', 'S:S373P', 'S:S375F',"S:D405N","S:R408S"),
                # ('S:N969K','S:Q954H'),
                # ("S:S371F","S:S373P","S:S375F"),
                # ("S:K356T","S:S371F","S:S373P","S:S375F","S:T376A","S:L390F","S:R403K"),
                # ("S:K417N","S:N440K","S:V445P"),
                # ("S:H655Y","S:N679K","S:P681H","S:A653T"),
                # ("S:S477N","S:T478K","S:DEL483/483","S:E484K","S:F486P","S:Q498R","S:N501Y")
            ]
    
    for cluster0 in clusters:
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
        if l0<=1:
            print(f'No detected stepwise evolution for {cluster0}')
            continue
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
        fn0 = ','.join(cluster0).replace('/','_')
        fig.tight_layout()
        plt.savefig(f'../descent_plots/ww_evo_seq{fn0}.pdf',transparent=True)
        plt.close('all')

if __name__ == "__main__":
    main()