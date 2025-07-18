""" Detect cryptic mutation clusters in wastewater sequencing data """

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import argparse
import sys, os
import pandas as pd
from tqdm import tqdm

from outbreak_data import outbreak_data as od
from outbreak_data import authenticate_user
from outbreak_tools import crumbs

# Window for clinical data query
START_DATE = "2020-01-01"
END_DATE = "2025-12-31"

FREYJA_BARCODES = "../sars2_metadata/usher_barcodes.feather"

parser = argparse.ArgumentParser(
    description="Detect cryptic mutation clusters in wastewater sequencing data"
)
parser.add_argument(
    "--covariants_dir", help="Directory containing coVar (linked mutations) output", type=str
)
parser.add_argument(
    "--output", help="Output file name", default="cryptic_variants.tsv"
)
parser.add_argument(
    "--max_clinical_count", default=10,
    help='Maximum number of clinical hits to consider a cluster "cryptic"',
)


# Hide print statements from API calls
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def parse_aa_muts(muts):
    """Prepare query for aa mutations"""
    output = []
    for m in muts.split(' '):
        if m == "Unknown": # Frameshift indels, most likely due to sequencing errors
            return []
        if m.split(":")[1][0] == m.split(":")[1][-1]: # Omit synonymous mutations
            continue
        if ":INS" in m: # Omit insertions
            continue
        if ":DEL" in m:
            if '/' not in m: # Workaround for single aa deletion query bug (e.g. S:DEL144 -> S:DEL144/144)
                output.append(f'{m}/{m.split(":DEL")[1]}')
            else:
                output.append(m)
        else:
            output.append(m)
    return list(set(output))


def parse_covariants(covariants_dir):
    """Parse covar output files, aggregate into one dataframe"""

    agg_covariants = pd.DataFrame()
    for file in os.listdir(covariants_dir):
        df = pd.read_csv(f'{covariants_dir}/{file}', sep="\t")

        df["query"] = df["aa_mutations"].apply(parse_aa_muts)
        df = df[df["query"].apply(len) > 0]
        df["sample"] = file
        agg_covariants = pd.concat([agg_covariants, df])

    return agg_covariants

def add_metadata(aggregate_covariants, metadata_file="../search_metadata.csv"):
    """Add metadata to aggregate covariants dataframe"""

    metadata = pd.read_csv(metadata_file)
    aggregate_covariants['sample'] = aggregate_covariants['sample'].str.split('.trimmed').str[0]

    # Get collection date and location from metadata
    aggregate_covariants = aggregate_covariants.merge(
        metadata[["sample", "collection_date", "location"]],
        on="sample",
        how="left",
    )

    aggregate_covariants["collection_date"] = pd.to_datetime(aggregate_covariants["collection_date"])

    # Drop rows with missing collection date or location
    aggregate_covariants = aggregate_covariants.dropna(subset=["collection_date", "location"])

    return aggregate_covariants


def query_clinical_data(aggregate_covariants, freyja_barcodes, START_DATE, END_DATE):
    """Query outbreak.info API for clinical detections of mutation clusters"""

    lineage_key = crumbs.get_alias_key()
    barcode_muts = pd.read_feather(freyja_barcodes).columns

    cache = {}
    for row in tqdm(aggregate_covariants.iterrows(), desc="Querying clinical data"):
        cluster = row[1]["query"]
        if str(cluster) in cache:
            continue
        if all([m in barcode_muts for m in row[1]["nt_mutations"].split(" ")]): # Skip if all mutations are in freyja barcodes
            cache[str(cluster)] = None
            continue

        try:
            with HiddenPrints():
                mut_data = od.lineage_cl_prevalence(
                    ".",
                    descendants=True,
                    mutations=cluster,
                    location="USA",
                    datemin=START_DATE,
                    datemax=END_DATE,
                    lineage_key=lineage_key,
                )
        except NameError as e:
            print(f"Error querying outbreak.info for cluster {cluster}: {e}")
            continue

        if mut_data is not None:
            cache[str(cluster)] = mut_data["lineage_count"].sum()
        else:
            cache[str(cluster)] = 0

    aggregate_covariants["num_clinical_detections"] = aggregate_covariants["query"].apply(
        lambda x: cache[str(x)] if str(x) in cache else None
    )

    aggregate_covariants = aggregate_covariants.dropna(subset=["num_clinical_detections"])

    return aggregate_covariants


def main():
    args = parser.parse_args()

    # Authenticate with GISAID credentials
    try:
        authenticate_user.get_authentication()
    except:
        authenticate_user.authenticate_new_user()

    aggregate_covariants = parse_covariants(args.covariants_dir)

    # Add metadata to aggregate covariants
    aggregate_covariants = add_metadata(aggregate_covariants)

    # Query clinical data
    cryptic_variants = query_clinical_data(
        aggregate_covariants, FREYJA_BARCODES, START_DATE, END_DATE
    )

    # Filter for cryptic variants
    cryptic_variants = cryptic_variants[
        cryptic_variants["num_clinical_detections"] <= float(args.max_clinical_count)
    ]

    # Save output
    cryptic_variants.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()