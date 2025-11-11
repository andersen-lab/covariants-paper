import os
import sys
import json
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from datetime import timedelta
from collections import Counter, defaultdict

RESULTS_DIR = "results"
FIGURE_DIR = "figures"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

genes = [
    {"name": "ORF1a", "start": 266, "end": 13483},
    {"name": "ORF1b", "start": 13468, "end": 21555},
    {"name": "S1",   "start": 21563, "end": 23926},
    {"name": "RBD",  "start": 22517, "end": 23185},
    {"name": "S2",   "start": 23927, "end": 25384},
    {"name": "E",    "start": 26245, "end": 26472},
    {"name": "M",    "start": 26523, "end": 27191},
    {"name": "N",    "start": 28274, "end": 29533},
]

def read_fasta_sequences(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    return {seq.id: str(seq.seq) for seq in sequences}


def load_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, parse_dates=["collection_date"])
    return metadata.sort_values("collection_date")


def generate_two_week_windows(metadata):
    two_week_windows, dates = [], []
    if metadata.empty:
        return two_week_windows, dates

    current_start = metadata["collection_date"].min()
    end_date = metadata["collection_date"].max()

    while current_start <= end_date:
        current_end = current_start + timedelta(weeks=2)
        mask = (metadata["collection_date"] >= current_start) & (
            metadata["collection_date"] < current_end
        )
        window_samples = metadata.loc[mask, "fasta_hdr"].tolist()
        if window_samples:
            two_week_windows.append(window_samples)
            dates.append((current_start, current_end))
        current_start = current_end

    return two_week_windows, dates


def load_fasta_to_srr(json_path):
    with open(json_path) as f:
        return json.load(f)


def convert_fasta_to_srr(two_week_windows, fasta_to_srr):
    srr_window, failures = [], 0
    for sample_window in two_week_windows:
        temp = []
        for sample in sample_window:
            if sample in fasta_to_srr:
                temp.append(fasta_to_srr[sample][0])
            else:
                failures += 1
        srr_window.append(temp)
    print("Number of missing mappings:", failures)
    return srr_window


def process_clinical_batch(sequences):
    position_freqs = defaultdict(dict)
    for i in range(len(sequences[0])):
        if i < 300 or i > 28800:
            continue
        column = [seq[i] for seq in sequences if seq[i].lower() != "n"]
        counts = Counter(column)
        total = sum(counts.values())
        position_freqs[i] = {b: c / total for b, c in counts.items()}
    return position_freqs


def entropy_per_position(position_freqs):
    def shannon(p):
        return -sum(v * math.log2(v) for v in p.values() if v > 0)

    return {pos: shannon(freqs) for pos, freqs in position_freqs.items()}

def aggregate_entropy(entropy_dict):
    return sum(entropy_dict.values())

def calculate_frequencies_window(srr_window, clinical_seqs, date_windows, output_path="./results/clinical_window_frequencies.json"):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Map keys like "ABC_12345" → "12345"
    processed = {key.split("_")[1]: seq for key, seq in clinical_seqs.items()}

    results_for_json = {}
    entropies = []

    for i, dataset in enumerate(srr_window):
        # Convert to a label like "2023-01-01_to_2023-01-15"
        start, end = date_windows[i]
        if not isinstance(start, str):
            start = start.strftime("%Y-%m-%d")
        if not isinstance(end, str):
            end = end.strftime("%Y-%m-%d")
        window_label = f"{start}_to_{end}"

        window_seqs = [processed[srr] for srr in dataset if srr in processed]

        if not window_seqs:
            results_for_json[window_label] = {}
            entropies.append(0)
            continue

        # Step 1: compute per-position frequencies
        position_freqs = process_clinical_batch(window_seqs)

        # Step 2: per-position entropy
        pos_entropies = entropy_per_position(position_freqs)

        # Step 3: aggregate entropy per window
        agg_entropy = aggregate_entropy(pos_entropies)
        entropies.append(agg_entropy)

        # Step 4: Save structure per window
        window_dict = {}
        for pos, freqs in position_freqs.items():
            pos_str = str(pos)
            window_dict[pos_str] = {**freqs}  # copy frequencies
            window_dict[pos_str]["Entropy"] = pos_entropies.get(pos, None)

        results_for_json[window_label] = window_dict

    # Save JSON file
    with open(output_path, "w") as f:
        json.dump(results_for_json, f, indent=4)

    return entropies


def parse_variants(filename):
    usecols = ['POS', 'REF', 'ALT', 'ALT_QUAL', 'ALT_FREQ', 'REF_DP', 'TOTAL_DP']
    df = pd.read_table(filename, usecols=usecols)
    df = df[
        (df['POS'] >= 300) & (df['POS'] <= 28800) &
        (df['ALT_QUAL'] >= 20) & (df['TOTAL_DP'] >= 10) &
        (df['ALT_FREQ'] > 0.001)
    ]
    df['REF_FREQ'] = df['REF_DP'] / df['TOTAL_DP']

    pos_dict = {}
    for row in df.itertuples(index=False):
        pos = row.POS
        alt = row.ALT
        alt_freq = row.ALT_FREQ
        ref = row.REF
        ref_freq = row.REF_FREQ
        if pos not in pos_dict:
            pos_dict[pos] = {}
            pos_dict[pos][ref] = ref_freq
        pos_dict[pos][alt] = alt_freq
    return pos_dict


def entropy_ww_samples(df, variant_path, suffix=".trimmed.sorted.unfiltered.tsv"):
    entropies, locations, dates, spike_entropies, non_spike_entropies = [], [], [], [], []
    site_specific_entropies = {}  # per sample
    for _, row in df.iterrows():
        sample = row['sample']
        filename = os.path.join(variant_path, sample + suffix)
        if not os.path.isfile(filename):
            filename = os.path.join(variant_path, sample + ".trimmed.sorted.tsv")
            if not os.path.isfile(filename):
                continue
        pos_freqs = parse_variants(filename)
        site_specific_entropies[sample] = pos_freqs

        e_per_position = entropy_per_position(pos_freqs)
        spike_dict = {}
        non_spike_dict = {}
        for key, value in e_per_position.items():
            if 21563 <= key <= 25384:
                spike_dict[key] = value
            else:
                non_spike_dict[key] = value

        spike_entropies.append(aggregate_entropy(spike_dict))
        non_spike_entropies.append(aggregate_entropy(non_spike_dict))
        entropies.append(aggregate_entropy(e_per_position))
        locations.append(row['location'])
        dates.append(row['collection_date'])
    return entropies, locations, dates, spike_entropies, non_spike_entropies, site_specific_entropies


def count_intrahost(srr_window, clinical_variants_path):
    counts_per_window = []
    positions_per_window = []
    for window in srr_window:
        window_positions = set()
        window_count = 0
        for filename in window:
            full_path = os.path.join(clinical_variants_path, filename + ".tsv")
            try:
                df = pd.read_table(
                    full_path,
                    usecols=["TOTAL_DP", "ALT_QUAL", "ALT_FREQ", "POS", "ALT"],
                    dtype={
                        "TOTAL_DP": "Int64",
                        "ALT_QUAL": "float64",
                        "ALT_FREQ": "float64",
                        "POS": "int32",
                        "ALT": "string"
                    }
                )
            except Exception:
                continue

            filtered_df = df[
                (df["TOTAL_DP"] >= 10) &
                (df["ALT_QUAL"] >= 20) &
                (df["ALT_FREQ"] >= 0.05) &
                (df["ALT_FREQ"] < 0.50) &
                (df["POS"] > 300) &
                (df["POS"] < 28800) &
                (~df["ALT"].str.contains(r"[\+\-]", na=False))
            ]
            window_count += len(np.unique(filtered_df['POS']))
            window_positions.update(filtered_df['POS'].unique())

        counts_per_window.append(window_count)
        positions_per_window.append(sorted(list(window_positions)))
    return counts_per_window, positions_per_window


def plot_clinical_vs_wastewater(clinical_df, wastewater_df, ihv_df, output_file=None):
    sns.set(style="whitegrid", font_scale=1.1)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    sns.lineplot(data=clinical_df, x="Midpoint", y="Entropy", color="tab:blue", marker="o", ax=ax1, label="Entropy")
    ax1b = ax1.twinx()
    sns.lineplot(data=ihv_df, x="Midpoint", y="Mean_IHV_Count", color="tab:red", linestyle="--", marker="s", ax=ax1b, label="IHV Count")

    ax1.set_ylabel("Shannon Entropy", color="tab:blue")
    ax1b.set_ylabel("IHV Count", color="tab:red")
    ax1.legend(loc="upper left")
    ax1b.legend(loc="upper right")

    sns.lineplot(data=wastewater_df, x="Date", y="Entropy", hue="Location", marker="o", ax=ax2)
    ax2.set_xlabel("Date")
    ax2.set_ylabel("Shannon Entropy")
    ax2.legend(loc="upper left", title="")

    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=300)
    plt.close()


def plot_entropy_by_protein(ww_df, output_file=None):
    sns.set_style("whitegrid")
    unique_locations = ww_df['Location'].unique()
    n_locations = len(unique_locations)

    fig, axes = plt.subplots(n_locations, 1, figsize=(12, 3 * n_locations), sharex=True)
    if n_locations == 1:
        axes = [axes]

    for ax, loc in zip(axes, unique_locations):
        subset = ww_df[ww_df['Location'] == loc]
        sns.lineplot(data=subset, x="Date", y="Shannon Entropy", hue="Type", ax=ax, linewidth=1.5)
        ax.set_title(f"Location: {loc}", loc='left')
        ax.set_ylabel("Shannon Entropy")
        ax.legend(title="Type", loc='upper right')

    axes[-1].set_xlabel("Date")
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=300)
    plt.close()

def calculate_and_save_all_results():
    # Paths
    metadata_path = "filtered_metadata.csv"
    ww_metadata_path = "search_metadata.csv"
    fasta_to_srr_path = "fasta_to_srr.json"
    clinical_multifasta = "clinical_multifasta.aligned.fa"
    ww_variant_path = "./wastewater_variants"
    clinical_variants_path = "./var_files"

    metadata = load_metadata(metadata_path)
    ww_metadata = pd.read_csv(ww_metadata_path)

    """
    # Wastewater entropies
    ww_entropies, ww_locations, ww_dates, spike_entropies, non_spike_entropies, ww_site_specific = entropy_ww_samples(ww_metadata, ww_variant_path)

    pd.DataFrame({
        "Date": ww_dates,
        "Location": ww_locations,
        "Entropy": ww_entropies,
        "Spike_Entropy": spike_entropies,
        "Non_Spike_Entropy": non_spike_entropies
    }).to_csv(os.path.join(RESULTS_DIR, "wastewater_entropy.csv"), index=False)

    # Save site-specific entropies as JSON
    json_safe = {k: {str(pos): val for pos, val in v.items()} for k, v in ww_site_specific.items()}
    with open(os.path.join(RESULTS_DIR, "wastewater_site_specific_entropies.json"), "w") as f:
        json.dump(json_safe, f, indent=2)
    """
    #clinical entropies 
    two_week_windows, dates = generate_two_week_windows(metadata)
    fasta_to_srr = load_fasta_to_srr(fasta_to_srr_path)
    srr_window = convert_fasta_to_srr(two_week_windows, fasta_to_srr)
    #clinical_seqs = read_fasta_sequences(clinical_multifasta)


    #aggregate_wastewater_entropy(ww_metadata, dates, ww_variant_path)
    #clinical_entropies = calculate_frequencies_window(srr_window, clinical_seqs, dates)
    """
    pd.DataFrame({
        "Start_Date": [d[0] for d in dates],
        "End_Date": [d[1] for d in dates],
        "Entropy": clinical_entropies
    }).to_csv(os.path.join(RESULTS_DIR, "clinical_entropy.csv"), index=False)
    """
    # Intrahost variant counts
    intrahost_counts, intrahost_positions = count_intrahost(srr_window, clinical_variants_path)
    pd.DataFrame({
        "Start_Date": [d[0] for d in dates],
        "End_Date": [d[1] for d in dates],
        "Mean_IHV_Count": intrahost_counts,
        "IHV_Positions": [",".join(map(str, pos_list)) for pos_list in intrahost_positions]
    }).to_csv(os.path.join(RESULTS_DIR, "intrahost_variant_counts.csv"), index=False)


def plot_from_saved_results():
    # Load
    ww_df = pd.read_csv(os.path.join(RESULTS_DIR, "wastewater_entropy.csv"))
    clinical_df = pd.read_csv(os.path.join(RESULTS_DIR, "clinical_entropy.csv"))
    ihv_df = pd.read_csv(os.path.join(RESULTS_DIR, "intrahost_variant_counts.csv"))

    # Midpoints for plotting
    clinical_df['Midpoint'] = pd.to_datetime(clinical_df['Start_Date']) + \
                              (pd.to_datetime(clinical_df['End_Date']) - pd.to_datetime(clinical_df['Start_Date'])) / 2
    ihv_df['Midpoint'] = pd.to_datetime(ihv_df['Start_Date']) + \
                         (pd.to_datetime(ihv_df['End_Date']) - pd.to_datetime(ihv_df['Start_Date'])) / 2
    ww_df['Date'] = pd.to_datetime(ww_df['Date'])

    # Filter by cutoff
    cutoff = pd.Timestamp("2024-02-01")
    clinical_plot_df = clinical_df[clinical_df['Midpoint'] < cutoff]
    wastewater_plot_df = ww_df[ww_df['Date'] < cutoff]
    ihv_plot_df = ihv_df[ihv_df['Midpoint'] < cutoff]

    # Call plotting functions
    plot_clinical_vs_wastewater(clinical_plot_df, wastewater_plot_df, ihv_plot_df,
                                output_file=os.path.join(RESULTS_DIR, "clinical_vs_wastewater_entropy.png"))


    ww_long = ww_df.melt(id_vars=["Date", "Location"], 
                         value_vars=["Entropy", "Spike_Entropy", "Non_Spike_Entropy"],
                         var_name="Type", value_name="Shannon Entropy")

    print("plotting entropy per protein")                        
    plot_entropy_by_protein(ww_long, output_file=os.path.join(RESULTS_DIR, "entropy_by_protein.png"))

    print("plotting average site entropy")
    plot_average_site_entropy()

def aggregate_wastewater_entropy(ww_metadata, dates, variant_path, output_json="./results/agg_wastewater_entropy.json", suffixes=(".trimmed.sorted.unfiltered.tsv", ".trimmed.sorted.tsv")): 
    results = {}  
    if not np.issubdtype(ww_metadata["collection_date"].dtype, np.datetime64): 
        ww_metadata["collection_date"] = pd.to_datetime(ww_metadata["collection_date"]) 
    for start, end in dates: 
        mask = (ww_metadata["collection_date"] >= start) & (ww_metadata["collection_date"] < end) & (ww_metadata["location"] == "Point Loma")
        window_samples = ww_metadata.loc[mask, "sample"].tolist() 

        if not window_samples: 
            continue 
        for sample in window_samples: 
            file_found = False 
            for suffix in suffixes: 
                filename = os.path.join(variant_path, sample + suffix) 
                if os.path.isfile(filename): 
                    file_found = True 
                    break 
            if not file_found: 
                continue 
        df = pd.read_table(filename, usecols=['POS', 'REF', 'ALT', 'ALT_DP', 'TOTAL_DP', 'ALT_QUAL', 'REF_DP']) 
        df = df[(df['POS'] >= 300) & (df['POS'] <= 28800) & (df['ALT_QUAL'] >= 20)] 

        agg_counts = {}
        for _, row in df.iterrows(): 
            pos, alt, ref = int(row['POS']), row['ALT'], row["REF"]
            if "+" in alt:
                continue
            if "-" in alt:
                continue 
            if pos not in agg_counts:
                agg_counts[pos] = {"TOTAL_DP":0}
            if alt not in agg_counts[pos]:
                agg_counts[pos][alt] = {'DP': 0}
            if ref not in agg_counts[pos]:
                agg_counts[pos][ref] = {'DP': 0 }
            
            agg_counts[pos][alt]['DP'] += row['ALT_DP'] 
            agg_counts[pos]['TOTAL_DP'] += row['TOTAL_DP'] 
            agg_counts[pos][ref]['DP'] += row['REF_DP']

        position_freqs = {} 
        for pos, var_dict in agg_counts.items(): 
            freqs = {}
            total_reads = var_dict['TOTAL_DP'] 
            if total_reads == 0: 
                continue 
            for k, v in var_dict.items():
                if k == "TOTAL_DP":
                    continue
                freqs[k] = v['DP'] / total_reads
            position_freqs[pos] = freqs

        pos_entropy = entropy_per_position(position_freqs) 
        window_key = f"{start.date()}_to_{end.date()}"
        results[window_key] = {
            str(pos): {**freqs, "Entropy": pos_entropy.get(pos, None)}
            for pos, freqs in position_freqs.items()
        }
    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(results, f, indent=2)

def normalize_frequencies(freq_dict):
    freq_dict = {k: v for k, v in freq_dict.items() if '+' not in k and "-" not in k}
    total = sum(freq_dict.values())
    if total == 0:
        return freq_dict  # avoid division by zero
    return {k: v / total for k, v in freq_dict.items()}

def plot_average_site_entropy(json_file=os.path.join(RESULTS_DIR, "wastewater_site_specific_entropies.json"),
                              output_file="average_gene_entropy_heatmap.png"):

    # Load data
    with open(json_file) as f:
        site_data = json.load(f)

    # Flatten JSON to DataFrame
    records = []
    for sample, pos_dict in site_data.items():
        #if sample != "1.05.23.PLJAN01.R1__NA__NA__230107_WW__00X":
        #    continue
        for pos_str, freqs in pos_dict.items():
            records.append({'Sample': sample, 'Position': int(pos_str), 'Freqs': freqs})

    df = pd.DataFrame(records)

    # Shannon entropy per site
    def row_entropy(freq_dict):
        norm_freqs = normalize_frequencies(freq_dict)
        return -sum(v * math.log2(v) for v in norm_freqs.values() if v > 0)

    df['Entropy'] = df['Freqs'].apply(row_entropy)

    # Set seaborn style
    sns.set(style="whitegrid")
    gene_entropy = []
    for gene in genes:
        mask = (df['Position'] >= gene['start']) & (df['Position'] <= gene['end'])
        mean_entropy = df.loc[mask, 'Entropy'].mean()
        gene_entropy.append({'Gene': gene['name'], 'Mean_Entropy': mean_entropy})

    gene_entropy_df = pd.DataFrame(gene_entropy)

    # --- Plot gene-wise mean entropy as a barplot ---
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 5))

    sns.barplot(
        data=gene_entropy_df,
        x='Gene',
        y='Mean_Entropy',
        color='orange',
        edgecolor='black',
        ax=ax
    )

    ax.set_title("Mean Entropy Per Gene")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Mean Entropy")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

def plot_entropy_per_gene(
    wastewater_file="./results/agg_wastewater_entropy.json",
    clinical_file="./results/clinical_window_frequencies.json",
    output_file="entropy_per_gene.png"):


    def load_entropy_by_gene(json_file):
        """Returns a DataFrame: each row is (Gene, Entropy) across all time windows."""
        with open(json_file, 'r') as f:
            data = json.load(f)

        rows = []
        for window, positions in data.items():  # window = "YYYY-MM-DD_to_YYYY-MM-DD"
            for pos, values in positions.items():
                pos = int(pos)
                entropy = values.get("Entropy", None)
                if entropy is None:
                    continue
                for gene in genes:
                    if gene["start"] <= pos <= gene["end"]:
                        rows.append({"Gene": gene["name"], "Entropy": entropy})
        return pd.DataFrame(rows)

    # Load
    ww_df = load_entropy_by_gene(wastewater_file)
    clinical_df = load_entropy_by_gene(clinical_file)

    # Average entropy per gene
    ww_avg = ww_df.groupby("Gene")["Entropy"].mean().reset_index()
    clinical_avg = clinical_df.groupby("Gene")["Entropy"].mean().reset_index()

    # Order genes on x-axis as defined
    gene_order = [g["name"] for g in genes]

    # Plotting
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)

    # Wastewater entropy
    sns.barplot(
        data=ww_avg,
        x="Gene",
        y="Entropy",
        order=gene_order,
        ax=axes[0],
        color="steelblue"
    )
    axes[0].set_title("Average Entropy per Gene — Wastewater Samples")
    axes[0].set_ylabel("Mean Entropy")

    # Clinical entropy
    sns.barplot(
        data=clinical_avg,
        x="Gene",
        y="Entropy",
        order=gene_order,
        ax=axes[1],
        color="salmon"
    )
    axes[1].set_title("Average Entropy per Gene — Clinical Samples")
    axes[1].set_xlabel("Genomic Region")
    axes[1].set_ylabel("Mean Entropy")

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved figure as {output_file}")
    
def plot_entropy_per_site(
    wastewater_file="./results/agg_wastewater_entropy.json", clinical_file="./results/clinical_window_frequencies.json", output_file="entropy_per_site.png"):

    def mean_entropy_per_pos(json_file):
        """Return DataFrame with Position, Mean_Entropy across all time windows."""
        with open(json_file, 'r') as f:
            data = json.load(f)

        rows = []
        for window, pos_data in data.items():  # window: "YYYY-MM-DD_to_YYYY-MM-DD"
            for pos, values in pos_data.items():
                entropy = values.get("Entropy", None)
                if entropy is not None:
                    rows.append({"Position": int(pos), "Entropy": entropy})

        df = pd.DataFrame(rows)
        return df.groupby("Position")["Entropy"].mean().reset_index().rename(columns={"Entropy": "Mean_Entropy"})

    # Load mean entropies from each dataset
    ww_entropy = mean_entropy_per_pos(wastewater_file).rename(columns={"Mean_Entropy": "Wastewater_Entropy"})
    clinical_entropy = mean_entropy_per_pos(clinical_file).rename(columns={"Mean_Entropy": "Clinical_Entropy"})

    # Merge on positions
    merged = pd.merge(ww_entropy, clinical_entropy, on="Position", how="inner")
    for index, row in merged.iterrows():
        ww = row['Wastewater_Entropy']
        cl = row['Clinical_Entropy']
        if abs(ww-cl) > 0.10:
            print(row['Position'], ww, cl)

    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=merged,
        x="Clinical_Entropy",
        y="Wastewater_Entropy",
        s=15,
        alpha=0.7
    )
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    lims = [max(min(xmin, ymin), 0), max(xmax, ymax)]  # Ensure square axes if needed

    plt.plot(lims, lims, '--', linewidth=1)
    plt.xlim(lims)
    plt.ylim(lims)
    plt.title("Per-Site Entropy: Clinical vs Wastewater")
    plt.xlabel("Clinical Mean Entropy")
    plt.ylabel("Wastewater Mean Entropy")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved scatterplot as {output_file}")

def plot_entropy_over_time_Sgene(
    wastewater_file="./results/agg_wastewater_entropy.json",
    clinical_file="./results/clinical_window_frequencies.json",
    output_file="entropy_over_time_Sgene.png"
):
    """
    Plot wastewater and clinical mean entropy over time using dual y-axes,
    only considering positions within the S gene (21563-25384).
    """

    S_START = 21563
    S_END = 25384

    def mean_entropy_per_window_Sgene(json_file):
        """Return DataFrame with window midpoint and mean entropy across S-gene positions only."""
        with open(json_file, 'r') as f:
            data = json.load(f)

        rows = []
        for window_label, positions in data.items():  # e.g., "2023-01-01_to_2023-01-15"
            start_str, end_str = window_label.split("_to_")
            start = pd.to_datetime(start_str)
            end = pd.to_datetime(end_str)
            midpoint = start + (end - start) / 2

            # Filter positions to S gene
            s_positions = [v["Entropy"] for pos, v in positions.items() 
                           if "Entropy" in v and S_START <= int(pos) <= S_END]

            if s_positions:
                mean_entropy = sum(s_positions) / len(s_positions)
                rows.append({"Time": midpoint, "Mean_Entropy": mean_entropy})

        df = pd.DataFrame(rows).sort_values("Time")
        return df

    # Load filtered data
    ww_df = mean_entropy_per_window_Sgene(wastewater_file).rename(columns={"Mean_Entropy": "Wastewater"})
    clinical_df = mean_entropy_per_window_Sgene(clinical_file).rename(columns={"Mean_Entropy": "Clinical"})

    # Plot
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Wastewater on left y-axis
    ax1.plot(ww_df["Time"], ww_df["Wastewater"], color="steelblue", marker="o", label="Wastewater")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Wastewater Mean Entropy", color="steelblue")
    ax1.tick_params(axis='y', labelcolor="steelblue")

    # Clinical on right y-axis
    ax2 = ax1.twinx()
    ax2.plot(clinical_df["Time"], clinical_df["Clinical"], color="salmon", marker="s", label="Clinical")
    ax2.set_ylabel("Clinical Mean Entropy", color="salmon")
    ax2.tick_params(axis='y', labelcolor="salmon")

    # Title and formatting
    plt.title("Wastewater vs Clinical Mean Entropy Over Time (S Gene)")
    fig.autofmt_xdate(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved figure as {output_file}")

def plot_normalized_entropy_per_gene(
    wastewater_file="./results/agg_wastewater_entropy.json",
    clinical_file="./results/clinical_window_frequencies.json",
    intrahost_file=os.path.join("results", "intrahost_variant_counts.csv"),
    output_file="normalized_entropy_per_gene.png",
    FIGURE_DIR="./figures"):

    def load_entropy_by_gene(json_file):
        """Returns a DataFrame: each row is (Gene, Entropy) across all time windows."""
        with open(json_file, 'r') as f:
            data = json.load(f)

        rows = []
        for window, positions in data.items():  # window = "YYYY-MM-DD_to_YYYY-MM-DD"
            for pos, values in positions.items():
                pos = int(pos)
                entropy = values.get("Entropy", None)
                if entropy is None:
                    continue
                for gene in genes:
                    if gene["start"] <= pos <= gene["end"]:
                        rows.append({"Gene": gene["name"], "Entropy": entropy})
        return pd.DataFrame(rows)

    # --- Load data ---
    ww_df = load_entropy_by_gene(wastewater_file)
    clinical_df = load_entropy_by_gene(clinical_file)
    ihv_df = pd.read_csv(intrahost_file)

    # --- Average entropy per gene ---
    ww_avg = ww_df.groupby("Gene")["Entropy"].mean().reset_index()
    clinical_avg = clinical_df.groupby("Gene")["Entropy"].mean().reset_index()

    # --- Normalize entropy values (min-max per dataset) ---
    for df in [ww_avg, clinical_avg]:
        min_val = df["Entropy"].min()
        max_val = df["Entropy"].max()
        df["Entropy_norm"] = (df["Entropy"] - min_val) / (max_val - min_val)

    # --- Process intrahost variant positions ---
    all_positions = []
    for pos_str in ihv_df["IHV_Positions"].dropna():
        pos_list = [int(p) for p in str(pos_str).split(",") if p.strip().isdigit()]
        all_positions.extend(pos_list)

    pos_counts = Counter(all_positions)

    # Map position counts to genes
    gene_counts = {gene["name"]: 0 for gene in genes}
    for pos, count in pos_counts.items():
        for gene in genes:
            if gene["start"] <= pos <= gene["end"]:
                gene_counts[gene["name"]] += count

    ihv_gene_data = []
    for gene in genes:
        gene_name = gene["name"]
        gene_length = gene["end"] - gene["start"] + 1
        raw_count = gene_counts.get(gene_name, 0)
        normalized_count = raw_count / gene_length
        ihv_gene_data.append({
            "Gene": gene_name,
            "IHV_Count": raw_count,
            "Gene_Length": gene_length,
            "IHV_per_base": normalized_count
        })

    ihv_gene_df = pd.DataFrame(ihv_gene_data)
    gene_order = [g["name"] for g in genes]

    # --- Plotting ---
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 12), sharex=True)

    # Wastewater normalized entropy
    sns.barplot(
        data=ww_avg,
        x="Gene",
        y="Entropy_norm",
        order=gene_order,
        ax=axes[0],
        color="steelblue"
    )
    axes[0].set_title("Normalized Entropy per Gene — Wastewater Samples")
    axes[0].set_ylabel("Normalized Entropy (0–1)")

    # Clinical normalized entropy
    sns.barplot(
        data=clinical_avg,
        x="Gene",
        y="Entropy_norm",
        order=gene_order,
        ax=axes[1],
        color="salmon"
    )
    axes[1].set_title("Normalized Entropy per Gene — Clinical Samples")
    axes[1].set_ylabel("Normalized Entropy (0–1)")

    # Intrahost variant counts normalized by gene length
    sns.barplot(
        data=ihv_gene_df,
        x="Gene",
        y="IHV_per_base",
        order=gene_order,
        ax=axes[2],
        color="mediumpurple"
    )
    axes[2].set_title("Intrahost Variant Density per Gene (Normalized by Gene Length)")
    axes[2].set_xlabel("Genomic Region")
    axes[2].set_ylabel("IHV per base")

    plt.xticks(rotation=45)
    plt.tight_layout()
    os.makedirs(FIGURE_DIR, exist_ok=True)
    plt.savefig(os.path.join(FIGURE_DIR, output_file), dpi=300)
    plt.close()
    print(f"Saved normalized figure as {os.path.join(FIGURE_DIR, output_file)}")


def plot_gene_dynamics_over_time(
    wastewater_file="./results/agg_wastewater_entropy.json",
    clinical_file="./results/clinical_window_frequencies.json",
    intrahost_file=os.path.join("results", "intrahost_variant_counts.csv"),
    output_file="gene_dynamics_over_time.png",
    FIGURE_DIR="./figures",
):
    """
    For each data modality (wastewater, clinical, intrahost),
    plot the per-gene metric across time:
      - Wastewater/clinical: mean entropy per gene per time window.
      - Intrahost: number of variants per gene normalized by gene length per window.
    """

    import json
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from collections import Counter
    import os

    # Helper to load entropy as (Window, Gene, Entropy)
    def load_entropy_by_gene_over_time(json_file):
        with open(json_file, "r") as f:
            data = json.load(f)
        rows = []
        for window, positions in data.items():  # window = "YYYY-MM-DD_to_YYYY-MM-DD"
            for pos, values in positions.items():
                pos = int(pos)
                entropy = values.get("Entropy", None)
                if entropy is None:
                    continue
                for gene in genes:
                    if gene["start"] <= pos <= gene["end"]:
                        rows.append({
                            "Window": window,
                            "Gene": gene["name"],
                            "Entropy": entropy
                        })
                        break
        return pd.DataFrame(rows)

    # --- Load datasets ---
    ww_df = load_entropy_by_gene_over_time(wastewater_file)
    clinical_df = load_entropy_by_gene_over_time(clinical_file)
    ihv_df = pd.read_csv(intrahost_file)

    # --- Aggregate mean entropy per gene per window ---
    ww_mean = ww_df.groupby(["Window", "Gene"])["Entropy"].mean().reset_index()
    clinical_mean = clinical_df.groupby(["Window", "Gene"])["Entropy"].mean().reset_index()

    # --- Process intrahost variant data ---
    # Expect columns: Start_Date, End_Date, Mean_IHV_Count, IHV_Positions
    rows = []
    for _, row in ihv_df.iterrows():
        window = f"{row['Start_Date']}_to_{row['End_Date']}"
        pos_list = []
        if pd.notna(row.get("IHV_Positions")):
            pos_list = [int(p) for p in str(row["IHV_Positions"]).split(",") if p.strip().isdigit()]

        pos_counts = Counter(pos_list)
        for pos, count in pos_counts.items():
            for gene in genes:
                if gene["start"] <= pos <= gene["end"]:
                    norm_count = count / (gene["end"] - gene["start"] + 1)
                    rows.append({
                        "Window": window,
                        "Gene": gene["name"],
                        "IHV_Count_norm": norm_count
                    })
                    break
    ihv_gene_df = pd.DataFrame(rows)

    # --- Aggregate intrahost counts per gene per window ---
    ihv_sum = ihv_gene_df.groupby(["Window", "Gene"])["IHV_Count_norm"].sum().reset_index()

    # --- Prepare plot ---
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 14), sharex=True)

    # Ensure window order is consistent (chronological)
    def sort_windows(df):
        df["Window_start"] = pd.to_datetime(df["Window"].str.split("_to_").str[0], errors="coerce")
        return df.sort_values("Window_start")

    ww_mean = sort_windows(ww_mean)
    clinical_mean = sort_windows(clinical_mean)
    ihv_sum = sort_windows(ihv_sum)
    window_order = sorted(ww_mean["Window"].unique(), key=lambda w: pd.to_datetime(w.split("_to_")[0]))

    # --- Wastewater ---
    sns.lineplot(
        data=ww_mean,
        x="Window",
        y="Entropy",
        hue="Gene",
        ax=axes[0],
        marker="o"
    )
    axes[0].set_title("Wastewater: Mean Entropy per Gene Over Time")
    axes[0].set_ylabel("Mean Entropy")
    axes[0].tick_params(axis='x', rotation=45)

    # --- Clinical ---
    sns.lineplot(
        data=clinical_mean,
        x="Window",
        y="Entropy",
        hue="Gene",
        ax=axes[1],
        marker="o"
    )
    axes[1].set_title("Clinical: Mean Entropy per Gene Over Time")
    axes[1].set_ylabel("Mean Entropy")
    axes[1].tick_params(axis='x', rotation=45)

    # --- Intrahost ---
    sns.lineplot(
        data=ihv_sum,
        x="Window",
        y="IHV_Count_norm",
        hue="Gene",
        ax=axes[2],
        marker="o"
    )
    axes[2].set_title("Intrahost: Normalized IHV Counts per Gene Over Time")
    axes[2].set_xlabel("Time Window")
    axes[2].set_ylabel("IHV Count per nt (normalized)")
    axes[2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    os.makedirs(FIGURE_DIR, exist_ok=True)
    plt.savefig(os.path.join(FIGURE_DIR, output_file), dpi=300)
    plt.close()

    print(f"Saved time-resolved figure as {os.path.join(FIGURE_DIR, output_file)}")

def main():
    calculate_and_save_all_results()

    plot_normalized_entropy_per_gene()
    plot_gene_dynamics_over_time()

    #plot_entropy_per_gene()
    #scatterplot of entropy where each point is a site and clinical vs wastewater
    #plot_entropy_per_site()

    #plot_entropy_over_time_Sgene()

    #print("Plotting from saved results...")
    #plot_from_saved_results()


if __name__ == "__main__":
    main()
