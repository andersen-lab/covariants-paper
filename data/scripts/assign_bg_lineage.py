import pandas as pd

def assign_lineage(nt_mutations, barcodes):
    mutations = nt_mutations.split(' ')
    mutations_filtered = [m for m in mutations if m in barcodes.columns]

    barcodes_filtered = barcodes[mutations_filtered + ['index']]
    barcodes_filtered.set_index('index', inplace=True)

    best_match = barcodes_filtered.sum(axis=1)

    return barcodes_filtered.loc[best_match].name
    

barcodes = pd.read_feather('../sars2_metadata/usher_barcodes.feather')
covar_output = pd.read_csv('../covar_clinical_detections.tsv', sep='\t')

covar_output['lineage'] = covar_output.apply(lambda row: assign_lineage(row['nt_mutations'], barcodes), axis=1)

covar_output.to_csv('../covar_clinical_detections_with_lineages.tsv', sep='\t', index=False)
