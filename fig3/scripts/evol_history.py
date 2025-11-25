import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib

from networkx.drawing.nx_agraph import graphviz_layout


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

MIN_LONGEST_PATH_LENGTH = 3
MAX_OUT_DEGREE = 1


def main():
    parser = argparse.ArgumentParser(description='Plot stepwise evolution of cryptic mutation clusters.')

    parser.add_argument('--input', type=str, help='Path to the input file (tsv) containing cryptic variants.')
    parser.add_argument('--output', type=str, default='../evoplots', help='Directory to save the output plots.')
    parser.add_argument('--min_muts', type=int, default=1, help='Minimum number of mutations in a cluster to consider.')
    parser.add_argument('--max_clinical_detections', type=int, default=100, help='Maximum number of clinical detections to consider a cluster as cryptic.')
    parser.add_argument('--min_observations', type=int, default=2, help='Minimum number of observations for a cluster to be included in the analysis.')
    parser.add_argument('--min_depth', type=int, default=10, help='Minimum depth for a mutation cluster to be considered.')

    args = parser.parse_args()

    df = pd.read_csv(args.input ,sep='\t')

    df['cluster'] = df['query'].apply(parse_query_list)
    df['cluster'] = df['cluster'].apply(lambda x: tuple(sorted(x,key=get_aa_site)))
    df['collection_date'] = pd.to_datetime(df['collection_date'])

    df['num_muts'] = df['cluster'].apply(lambda x:len(x))
    df = df[df['num_muts']>= args.min_muts]

    df.drop_duplicates(subset=['cluster', 'collection_date', 'location'],keep='first',inplace=True)

    df = df[df['cluster_depth'] >= args.min_depth]

    df['coverage_start'] = df['coverage_start'].apply(lambda x: (x - 21563) // 3)
    df['coverage_end'] = df['coverage_end'].apply(lambda x: (x - 21563) // 3)

    # aggregate metadata for each unique cluster
    df_aggregate = df.groupby('cluster').agg(
        {
            'cluster_depth':tuple,
            'location':tuple,
            'collection_date':tuple,
            'coverage_start':tuple,
            'coverage_end':tuple,
            'num_clinical_detections':'mean',
            'num_muts':'min'
        }
    ).sort_values(by='num_muts',ascending=False)
    
    df_aggregate['total_observations'] = df_aggregate['cluster_depth'].apply(lambda x: len(x))
    df_aggregate = df_aggregate[df_aggregate['total_observations'] >= args.min_observations]

    df_aggregate = df_aggregate[df_aggregate['num_clinical_detections'] <= args.max_clinical_detections]
    
    df_aggregate['cluster'] = df_aggregate.index

    test_clusters = []
    for row in df_aggregate.iterrows():
        test_clust = set(row[1]['cluster'])
        df_subset = df_aggregate[df_aggregate['cluster'].apply(lambda x: set(x).issubset(test_clust) & (set(x)!=test_clust))]
        test_clusters.append(df_aggregate.drop(index=df_subset.index)['cluster'].to_list())

    # flatten
    test_clusters = [item for sublist in test_clusters for item in sublist if len(item) > 0]
    test_clusters = list(set(test_clusters))

    df_aggregate.to_csv('df_aggregate.tsv', sep='\t')

    count = 0
    #generate subtrees leading up to these cryptics. 
    for test_clust in test_clusters:

        df_region = df_aggregate.copy()

        positions = np.array([get_aa_site(t) for t in test_clust])
        df_region = df_region[df_region['coverage_start'].apply(lambda x: all([pos <= min(positions) for pos in x]))]
        df_region = df_region[df_region['coverage_end'].apply(lambda x: all([pos >= max(positions) for pos in x]))]

        if len(df_region) == 0:
            continue

        df_superset = df_region[df_region['cluster'].apply(lambda x: set(x).issubset(test_clust))]
        df_superset.loc[:, 'cluster'] = df_superset['cluster'].apply(lambda x: tuple(sorted(list(x), key=get_aa_site)))
        parent_list = []
        new_muts_list = []
        for c in df_superset['cluster']:
            # get subsets
            subsets = df_superset[df_superset['cluster'].apply(lambda x: set(x).issubset(c) & (x != c))]
            if subsets.shape[0]>0:
                subsets = pd.concat((subsets,pd.Series(subsets['cluster'].apply(lambda x: len(set(c) - set(x))),name='difference')),axis=1)
                minDiff = subsets['difference'].min()
                subsets = subsets[subsets['difference'] == minDiff]
                
                if subsets.shape[0] == 1:
                    parent_list.append([subsets['cluster'].iloc[0]])
                    new_muts_list.append([tuple(set(c)-set(subsets['cluster'].iloc[0]))])
                else:
                    parent_list.append(list(subsets['cluster']))
                    new_muts_list.append([tuple(set(c) - set(t0)) for t0 in subsets['cluster']])
            else:
                parent_list.append(None)
                new_muts_list.append(c)

        # for each, define parent, unless smallest cluster. 
        edge_labels = {}
        G = nx.DiGraph()
        # assemble graph
        # make edges to first level
        for j,c in enumerate(df_superset['cluster']):
            if parent_list[j] is not None:
                for l in range(0,len(parent_list[j])):
                    G.add_edge(
                        parent_list[j][l],
                        c,
                        weight=1
                    )
                    edge_labels[(parent_list[j][l],c)] = ','.join(new_muts_list[j][l])        

        try:
            if len(nx.dag_longest_path(G)) <  MIN_LONGEST_PATH_LENGTH:
                continue
        except:
            pass

        fig, ax = plt.subplots()
        pos = graphviz_layout(
            G,
            prog = 'dot',
            args = "-Grankdir='LR' -Goverlap=false",
            root = 0
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            node_color = "cornflowerblue",
            alpha = 0.9,
            margins = 0.1,
            node_size = 80
        )

        nx.draw_networkx_edges(
            G,
            pos,
            width = 1,
            alpha = 0.8,
            arrows = True,
            edge_color = "grey",
        )

        labels = {c:o for c,o in zip(df_superset['cluster'],list(df_superset['total_observations']))}

        if len(labels) > len(pos):
            labels = {key:labels[key] for key in pos.keys()}
        
        nx.draw_networkx_labels(
            G,
            pos = pos,
            labels = labels,
            font_color = "white",
            font_size = 6
        )
        
        nx.draw_networkx_edge_labels(
            G,
            pos = pos,
            edge_labels = edge_labels,
            font_color = 'black',
            font_size = 5
        )

        ymin, ymax = ax.get_ylim()
        y_length = ymax - ymin
        
        xmin, xmax = ax.get_xlim()
        x_length = xmax - xmin

        # list mutations for starting and ending nodes
        for pk in pos.keys():
            if G.out_degree(pk) > MAX_OUT_DEGREE:
                continue
            if G.in_degree(pk)==0:
                ax.text(
                    pos[pk][0] - x_length * 0.06,
                    pos[pk][1],
                    '\n'.join(pk),
                    fontsize = 5,
                    verticalalignment = 'center',
                    horizontalalignment = 'right'
                )
            if G.out_degree(pk)==0:
                ax.text(
                    pos[pk][0] + x_length * 0.06,
                    pos[pk][1],
                    '\n'.join(pk),
                    fontsize = 5,
                    verticalalignment = 'center',
                    horizontalalignment ='left'
                )

        # add number of clinical detections below ww detection count
        for pk in pos.keys():
            ax.text(
                pos[pk][0],
                pos[pk][1] - y_length * 0.03,
                int(df_region.loc[[pk],'num_clinical_detections']),
                fontsize = 5,
                verticalalignment = 'center',
                horizontalalignment ='center',
                color = 'red'
            )
        # add the sites the variant was detected at. 
        for pk in pos.keys():
            locs = set(list(df_region.loc[[pk],'location'][0]))
            for j,l0 in enumerate(['South Bay','Point Loma','Encina']):
                if l0 in locs:
                    ax.scatter(
                        pos[pk][0] + x_length * 0.02,
                        pos[pk][1] + y_length * 0.02 * (1-j),
                        marker = 's',
                        s = 5,
                        color = 'darkblue',
                        linewidth = 0.5
                    )
                else:
                    ax.scatter(
                        pos[pk][0] + x_length * 0.02,
                        pos[pk][1] + y_length * 0.02 * (1-j),
                        marker = 's',
                        s = 5,
                        edgecolors = 'darkblue',
                        facecolors = 'none',
                        linewidth = 0.5
                    )
        
        for pk in pos.keys():
            tp = pd.to_datetime(pd.Series(df_region.loc[[pk],'collection_date'][0],name='times'))
            if len(tp)==1:
                    ax.text(
                        pos[pk][0],
                        pos[pk][1] + y_length * 0.05,
                        f"{tp.iloc[0].strftime('%Y-%m-%d')}",
                        fontsize = 4,
                        verticalalignment = 'center',
                        horizontalalignment = 'center',
                        color = 'black'
                    )
            else:
                ax.text(
                    pos[pk][0],
                    pos[pk][1] + y_length * 0.05,
                     f"{tp.min().strftime('%Y-%m-%d')}\n-{tp.max().strftime('%Y-%m-%d')}",
                    fontsize = 4,
                    verticalalignment = 'center',
                    horizontalalignment = 'center',
                    color = 'black'
                )

        plt.gca().set_frame_on(False)
        plt.savefig(f'{args.output}/ww_evo_seq{count}.pdf', dpi=300, transparent=True, bbox_inches='tight')
        count += 1
        plt.close('all')


def parse_query_list(query_list):
    """Parse the query list from the cryptic mutations file."""

    output = []
    query_list = query_list.strip('[]').split(', ')
    query_list = [x.strip().strip("'").strip('"') for x in query_list]
    for mut in query_list:
        if 'DEL' in mut:
            if mut.split('DEL')[1].split('/')[0] == mut.split('DEL')[1].split('/')[1]:
                output.append(mut.split('DEL')[0] + 'DEL' + mut.split('DEL')[1].split('/')[0])
            else:
                output.append(mut)
        else:
            output.append(mut)        
    return query_list

def get_aa_site(mut):
    """Extract the amino acid site from a mutation string."""

    if 'DEL' in mut:
        return int(mut.split('DEL')[1].split('/')[0])
    else: 
        return int(mut.split(':')[1][1:-1])

if __name__ == "__main__":
    main()