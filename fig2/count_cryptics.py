import pandas as pd

covar = pd.read_csv('../data/covar_clinical_detections.tsv', sep='\t')
covar = covar[covar['num_clinical_detections'] <= 10]

covar['collection_date'] = pd.to_datetime(covar['collection_date'])
covar['len_mutations'] = covar['nt_mutations'].apply(lambda x:len(x.split(' ')))
covar = covar[covar['len_mutations'] > 1]

# select only nt_mutations that occur more than once
covar_counts = covar.groupby('nt_mutations').size().reset_index(name='counts')
covar = covar.merge(covar_counts, on='nt_mutations', how='left')
covar = covar[covar['counts'] > 1]

# drop duplicates in nt_mutations
unique_cryptics = covar.drop_duplicates(subset=['nt_mutations'])

print(f"Total unique cryptic lineages: {len(unique_cryptics)}")


