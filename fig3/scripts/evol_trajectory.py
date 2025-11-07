import argparse
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

""" Given a cryptic mutation cluster, plot the stepwise evolution of the mutations that descent from it."""


def main():
    parser = argparse.ArgumentParser(description='Plot stepwise evolution of cryptic mutation clusters.')

    parser.add_argument('--input', type=str, help='Path to the input file (tsv) containing cryptic variants.')
    parser.add_argument('--output', type=str, default='../descent_plots', help='Directory to save the output plots.')
    parser.add_argument('--min_muts', type=int, default=2, help='Minimum number of mutations in a cluster to consider.')
    parser.add_argument('--max_clinical_detections', type=int, default=10, help='Maximum number of clinical detections to consider a cluster as cryptic.')
    parser.add_argument('--min_observations', type=int, default=2, help='Minimum number of observations for a cluster to be included in the analysis.')
    parser.add_argument('--min_freq', type=float, default=0.001, help='Minimum frequency of mutations in the cluster to be considered.')
    parser.add_argument('--min_depth', type=int, default=10, help='Minimum depth for a mutation cluster to be considered.')

    args = parser.parse_args()

    df = pd.read_csv(args.input ,sep='\t')

    df['cluster'] = df['query'].apply(parse_query_list)
    df['cluster'] = df['cluster'].apply(lambda x: tuple(sorted(x,key=get_aa_site)))
    df['collection_date'] = pd.to_datetime(df['collection_date'])

    df['num_muts'] = df['cluster'].apply(lambda x:len(x))
    df = df[df['num_muts'] >= args.min_muts]

    df.drop_duplicates(subset=['cluster', 'collection_date', 'location'],keep='first',inplace=True)

    df = df[df['cluster_depth'] >= args.min_depth]
    df = df[df['frequency'] >= args.min_freq]

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
    ).sort_values(by='num_muts', ascending=False)

    df_aggregate['total_observations']= df_aggregate['cluster_depth'].apply(lambda x: len(x))
    df_aggregate = df_aggregate[df_aggregate['total_observations'] >= args.min_observations]

    df_aggregate = df_aggregate[df_aggregate['num_clinical_detections'] <= args.max_clinical_detections]
    df_aggregate.to_csv('df_aggregate.tsv', sep='\t')
    df_aggregate['cluster'] = df_aggregate.index

    # for a known mutation cluster, see what evolves from it. 
    test_clusters = [
                ("S:K417N","S:N440K", "S:V445P", "S:G446S", "S:N460K"),#"S:V445P", "S:G446S", "S:N460K"), # XBB.1.5 (28.50%)
                # ("S:D614G", "S:H655Y"), # BA.1.1 (18.33%)
                # ('S:G142D', 'S:DEL144'), # XBB.1.5 (16.03%)
                # ("S:G142D", "S:DEL144", "S:F157S", "S:R158G"),
                # ##############################
                # ("S:S371F","S:S373P","S:S375F", "S:K356T", "S:T376A", "S:R403K"),
                # ("S:N764K"),
                # ("S:M697V", "S:N679K", "S:P681H"),
                # ("S:H655Y","S:T678A","S:N679K","S:P681R"),
                # ('S:K356T', 'S:S371F', 'S:S373P', 'S:S375F',"S:D405N","S:R408S"),
                # ('S:N969K','S:Q954H'),
                # ("S:K356T","S:S371F","S:S373P","S:S375F","S:T376A","S:L390F","S:R403K"),
                # ("S:H655Y","S:N679K","S:P681H","S:A653T"),
                # ("S:S477N","S:T478K")
            ]
    
    for test_clust in test_clusters:
        df_superset = df_aggregate[df_aggregate['cluster'].apply(lambda x: set(test_clust).issubset(set(x)) or set(test_clust)==set(x))]
        print(df_superset)

        df_superset['num_descendants'] = [df_superset[df_superset['cluster'].apply(lambda x: set(c0).issubset(set(x)))].shape[0] for c0 in df_superset.index] 
        df_superset = df_superset.sort_values(by='total_observations',ascending=False)

        # limit to 10 for simplicity
        if df_superset.shape[0] > 10:
            df_superset = df_superset.iloc[0:10]

        #force the seed cluster to be in the matrix
        if df_superset[df_aggregate['cluster'].apply(lambda x: set(test_clust)==set(x))].shape[0]==0:
            df_superset = df_superset.reset_index(drop=True)
            newRow = pd.DataFrame({c0:None for c0 in df_aggregate.columns},index=[test_clust])
            newRow['cluster'] =[test_clust]
            df_superset = pd.concat((df_superset,newRow),axis=0,ignore_index=True).set_index('cluster')
            df_superset['cluster'] = df_superset.index

        parent_list = []
        new_muts_list = []
        for j,c in enumerate(df_superset['cluster']):
            # get subsets
            subsets = df_superset[df_superset['cluster'].apply(lambda x: set(x).issubset(c) & (x!=c))]
            if subsets.shape[0]>0:
                subsets['difference'] = subsets['cluster'].apply(lambda x: len(set(c)-set(x)))
                minDiff = subsets['difference'].min()
                subsets = subsets[subsets['difference']==minDiff]
                
                if subsets.shape[0] ==1:
                    parent_list.append([subsets['cluster'].iloc[0]])
                    new_muts_list.append([tuple(set(c)-set(subsets['cluster'].iloc[0]))])
                else:
                    parent_list.append(list(subsets['cluster']))
                    new_muts_list.append([tuple(set(c)- set(t0)) for t0 in subsets['cluster']])
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
                    G.add_edge(parent_list[j][l],c,weight=1)
                    edge_labels[(parent_list[j][l],c)] = ','.join(new_muts_list[j][l])

        # Create time axis plot of the accumulated mutations using seaborn swarmplot
        plot_data = []
        for edge in edge_labels.keys():
            from_node, to_node = edge

            mut = edge_labels[edge]
            dates = pd.Series(df_superset.loc[[to_node], 'collection_date'].values[0])
            locations = pd.Series(df_superset.loc[[to_node], 'location'].values[0])
            
            for row in [{'Mutation': mut, 'Date': date, 'Location': loc} for date, loc in zip(dates, locations)]:
                plot_data.append(row)

        plot_df = pd.DataFrame.from_records(plot_data)
        plot_df['Date'] = pd.to_datetime(plot_df['Date'])

        custom_palette = {
            'South Bay': '#F7D842',
            'Point Loma': '#A8C763',
            'Encina': '#6FB8A6'
        }

        fig, ax = plt.subplots(figsize=(9, 11))
        fig = sns.swarmplot(
            data=plot_df,
            x='Date',
            y='Mutation',
            hue='Location',  # Color by Location
            ax=ax,
            size=5,
            palette=custom_palette
        )
        ax.legend(title='Location', bbox_to_anchor=(1.05, 1), loc='center left')

        fig = fig.get_figure()
        fig.savefig(f'{args.output}/mutation_swarmplot_{test_clust}.pdf', transparent=True, bbox_inches='tight')

        ax.set_xlabel('Collection Date')
        ax.set_ylabel('Mutation')
        ax.set_title('Cryptic Descendant Detections Over Time')
        plt.tight_layout()
        locator = mdates.MonthLocator(bymonthday=1)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.set_xlabel('Collection Date')
        ax.set_xlim(xmin=pd.Timestamp('2023-01-01'), xmax=pd.Timestamp('2024-01-01'))
        ax.spines[['right', 'top']].set_visible(False)


        plt.tight_layout()
        plt.savefig(f'{args.output}/mutation_swarmplot_{test_clust}.pdf', transparent=True, bbox_inches='tight')
        plt.close('all')

        l0 = len(nx.dag_longest_path(G))-1
        # if l0 <= 1:
        #     print(f'No detected stepwise evolution for {test_clust}')
        #     continue

        fig, ax = plt.subplots(figsize=(3*l0,5))

        pos = graphviz_layout(
            G,
            prog = 'dot',
            args="-Grankdir='LR' -Goverlap=false",
            root=0
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            node_color="cornflowerblue",
            alpha=0.9,
            margins=0.1,
            node_size=80
        )

        nx.draw_networkx_edges(
            G,
            pos,
            width=1,
            alpha=0.8,
            arrows=True,
            edge_color="grey",
        )

        labels = {c:o for c,o in zip(df_superset['cluster'],list(df_superset['total_observations']))}

        if len(labels) > len(pos):
            labels = {key:labels[key] for key in pos.keys()}

        labels.pop(test_clust)

        nx.draw_networkx_labels(
            G,
            pos=pos,
            labels=labels,
            font_color="white",
            font_size=6
        )

        nx.draw_networkx_edge_labels(
            G,
            pos=pos,
            edge_labels=edge_labels,
            font_color="black",
            font_size=5
        )

        ymin, ymax = ax.get_ylim()
        y_length = ymax - ymin

        xmin, xmax = ax.get_xlim()
        x_length = xmax - xmin

        for pk in pos.keys():
            if G.in_degree(pk) == 0:
                ax.text(
                    pos[pk][0] - x_length * 0.03,
                    pos[pk][1],
                    '\n'.join(pk),
                    fontsize=5,
                    verticalalignment='center',
                    horizontalalignment='right'
                )


        # add location markers
        for pk in pos.keys():
            if df_superset.loc[[pk],'location'][0] is not None:
                locs = set(list(df_superset.loc[[pk],'location'][0]))
            else:
                continue
            for j,l0 in enumerate(['South Bay','Point Loma','Encina']):
                if l0 in locs:
                    ax.scatter(
                        pos[pk][0] + x_length * 0.035,
                        pos[pk][1] + y_length * 0.015 * (1-j),
                        marker='s',
                        s=5,
                        color='darkblue',
                        linewidth=0.5
                    )
                else:
                    ax.scatter(
                        pos[pk][0] + x_length * 0.035,
                        pos[pk][1] + y_length * 0.015 * (1-j),
                        marker='s',
                        s=5,
                        edgecolors='darkblue',
                        facecolors='none',
                        linewidth=0.5
                    )
        
        # add number of clinical detections below ww detection count
        for pk in pos.keys():
            ax.text(
                pos[pk][0],
                pos[pk][1] - y_length * 0.03,
                str(int(df_superset.loc[[pk], 'num_clinical_detections'])) if pk != test_clust else "",
                fontsize=5,
                verticalalignment='center',
                horizontalalignment='center',
                color='red'
            )

        # Add collection date information
        for pk in pos.keys():
            if df_superset.loc[[pk], 'location'][0] is not None:
                tp = pd.Series(df_superset.loc[[pk], 'collection_date'][0], name='times')
            else:
                continue

            if len(tp) == 1:
                ax.text(
                    pos[pk][0],
                    pos[pk][1] + y_length * 0.03,
                    f"{tp.iloc[0].strftime('%Y/%m/%d')}",
                    fontsize=4,
                    verticalalignment='center',
                    horizontalalignment='center',
                    color='black'
                )
            else:
                ax.text(
                    pos[pk][0],
                    pos[pk][1] + y_length * 0.03,
                    f"{tp.min().strftime('%Y/%m/%d')}\n-{tp.max().strftime('%Y/%m/%d')}",
                    fontsize=4,
                    verticalalignment='center',
                    horizontalalignment='center',
                    color='black'
                )

        plt.gca().set_frame_on(False)
        fn0 = ','.join(test_clust).replace('/','_')
        fig.tight_layout()
        plt.savefig(f'{args.output}/ww_evo_seq{fn0}.pdf',transparent=True, bbox_inches='tight')
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
    return output


def get_aa_site(mut):
    if 'DEL' in mut:
        return int(mut.split('DEL')[1].split('/')[0])
    else: 
        return int(mut.split(':')[1][1:-1]
)


if __name__ == "__main__":
    main()