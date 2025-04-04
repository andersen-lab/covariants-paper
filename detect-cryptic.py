""" Detect cryptic mutation clusters in wastewater sequencing data """

import argparse
import os
import pandas as pd

from outbreak_data import outbreak_data as od
from outbreak_data import authenticate_user
from outbreak_tools import crumbs


def clean_raw_muts(muts):
    """Clean raw mutation strings from freyja output"""
    output = []
    for m in muts.split(" "):
        if ":" not in m or "*" in m:
            continue
        if ":INS" in m:
            continue
        if ":DEL" in m:
            deletion_size = m.split(")(")[0].split(",")[1]
            if int(deletion_size) % 3 != 0:
                continue
            output.append(m)
        else:
            output.append(m)
    if len(output) == 0:
        return []
    return output


def parse_aa_muts(muts):
    output = []
    for m in muts:
        if ":DEL" in m:
            if '/' not in m: # Workaround for deletion query bug
                output.append(f'{m.split(")(")[1][:-1]}/{m.split(":DEL")[1][:-1]}')
            else:
                output.append(m.split(")(")[1][:-1])
        else:
            output.append(m.split("(")[1][:-1])
    return list(set(output))


def parse_covariants(covariants_dir, metadata_file):
    """Parse freyja covariants output files, aggregate into one dataframe"""

    agg_covariants = pd.DataFrame()
    for file in os.listdir(covariants_dir):
        df = pd.read_csv(f'{covariants_dir}/{file}', sep="\t")
        try:
            df["Covariants"] = df["Covariants"].apply(clean_raw_muts)
        except:
            continue
        df["query"] = df["Covariants"].apply(parse_aa_muts)
        df = df[df["query"].apply(len) > 0]
        df["sample"] = file
        agg_covariants = pd.concat([agg_covariants, df])

    agg_covariants.to_csv("agg_covariants.tsv", sep="\t", index=False)

    # Merge metadata with covariants (if provided)
    if metadata_file is not None:
        metadata = pd.read_csv(metadata_file, sep="\t")
        agg_covariants = pd.merge(agg_covariants, metadata, on="sample", how="left")
    return agg_covariants


def query_clinical_data(aggregate_covariants, freyja_barcodes, START_DATE, END_DATE):
    """Query outbreak.info API for clinical detections of mutation clusters"""

    lineage_key = crumbs.get_alias_key()
    barcode_muts = pd.read_csv(freyja_barcodes, sep=",").columns
    cache = {}
    for row in aggregate_covariants.iterrows():
        cluster = row[1]["query"]
        if str(cluster) in cache:
            continue
        if all([m.split("(")[0] in barcode_muts for m in cluster]):
            cache[str(cluster)] = None
            continue

        try:
            mut_data = od.lineage_cl_prevalence(
                ".",
                descendants=True,
                mutations=cluster,
                location="USA",
                datemin=START_DATE,
                datemax=END_DATE,
                lineage_key=lineage_key,
            )
        except NameError:
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
    parser = argparse.ArgumentParser(
        description="Detect cryptic mutation clusters in wastewater sequencing data"
    )

    parser.add_argument(
        "--covariants_dir", help="Directory containing freyja covariants output", type=str
    )
    parser.add_argument(
        "--metadata",
        help='TSV containing at least "covariants_filename". Adds additional metadata to output file (see data/example_metadata.tsv)'
    )
    parser.add_argument(
        "--output", help="Output file name", default="cryptic_variants.tsv"
    )
    parser.add_argument(
        "--max_clinical_count",
        help='Maximum number of clinical samples to consider a cluster "cryptic"',
    )

    args = parser.parse_args()

    START_DATE = "2020-01-01"
    END_DATE = "2024-12-31"
    FREYJA_BARCODES = "freyja_metadata/sars-cov-2/usher_barcodes.csv"

    # Authenticate with GISAID credentials
    try:
        authenticate_user.get_authentication()
    except:
        authenticate_user.authenticate_new_user()

    aggregate_covariants = parse_covariants(args.covariants_dir, args.metadata)

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