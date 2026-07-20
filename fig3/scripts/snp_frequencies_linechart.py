"""Plot SNP frequencies over time as a line chart, binned by month."""
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIG3_DIR = SCRIPT_DIR.parent
DATA_DIR = FIG3_DIR.parent / "data"

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def main():
    freyja_barcodes = pd.read_feather(DATA_DIR / "sars2_metadata/usher_barcodes.feather")
    barcode_muts = freyja_barcodes.columns.tolist()

    df = pd.read_csv(DATA_DIR / "covar_clinical_detections.tsv", sep="\t")
    df["mutations_filtered"] = df["nt_mutations"].apply(
        lambda x: [mut for mut in x.split(" ") if mut not in barcode_muts]
    )

    # select only rows with mutations detected in multiple rows
    df["collection_date"] = pd.to_datetime(df["collection_date"])
    df["month"] = df["collection_date"].dt.to_period("M").dt.to_timestamp()

    # snp_types = [
    #     "A>T", "A>C", "A>G", "T>A", "T>C", "T>G",
    #     "C>A", "C>T", "C>G", "G>A", "G>C", "G>T",
    # ]

    snp_types = [
        "T>C", "A>G", "G>A", "C>T", "G>T"
    ]
    color_map = {"A>G": "#1f77b4", "T>C": "#ff7f0e", "C>T": "#2ca02c", "G>A": "#d62728", "G>T": "#9467bd"}
    colors = [color_map[st] for st in snp_types]

    # Compute SNP type counts per month
    monthly_snp_counts = []
    for _, row in df.iterrows():
        for mut in row["mutations_filtered"]:
            if "+" in mut or "-" in mut:
                continue
            snp_type = f"{mut[0]}>{mut[-1]}"
            if snp_type not in snp_types:
                continue
            monthly_snp_counts.append({"month": row["month"], "snp_type": snp_type})

    monthly_df = pd.DataFrame(monthly_snp_counts)
    monthly_freqs = monthly_df.groupby(["month", "snp_type"]).size().unstack(fill_value=0)
    #monthly_freqs = monthly_freqs.div(monthly_freqs.sum(axis=1), axis=0)


    for st in snp_types:
        if st not in monthly_freqs.columns:
            monthly_freqs[st] = 0
    monthly_freqs = monthly_freqs[snp_types]


    fig, ax = plt.subplots()
    for snp_type, color in zip(snp_types, colors):
        ax.plot(monthly_freqs.index, monthly_freqs[snp_type], label=snp_type, color=color)

    #ax.axhline(y=0.081, color="black", linestyle="--", alpha=0.7)
    #ax.set_ylim(0, 0.8)
    ax.set_title("Absolute Cryptic SNPs Counts Over Time (monthly Bins)")
    ax.set_xlabel("Collection Date (month)")
    ax.set_ylabel("Frequency")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", ncol=1)
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    out_path = FIG3_DIR / "cryptic_snp_absolutes_linechart.png"
    plt.savefig(out_path)
    plt.close()
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
