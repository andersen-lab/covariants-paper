# covariants-paper

Code and data for **"Read-level mutation linkage in wastewater sequencing reveals cryptic virus evolution and persistence"** 

The paper introduces **coVar** (https://github.com/andersen-lab/covar), a tool that identifies physically-linked mutations on sequencing read pairs, and combines it with the [outbreak.info](https://outbreak.info) API to flag mutation clusters that are rarely or never seen in global clinical sequences ("cryptic variants"). This repo contains the downstream analysis: turning coVar output from ~950 San Diego wastewater samples (Jan 2023–Dec 2024) into the paper's figures.

## How the pipeline fits together

Everything downstream starts from raw BAMs and ends at `data/covar_clinical_detections.tsv`, the single aggregated table that essentially every figure script reads:

```
data/search_bucket.txt
        │  download_search_data.sh
        ▼
   bam/*.bam ──────────────────────────────┐
        │  rename_reference.sh             │  get_median_coverage.sh (fig2/)
        │  (fixes contig naming)           │  avg_quality.sh (supp_fig3/)
        ▼                                  ▼
   parse_search_metadata.py    fig2/median_coverage_summary.tsv
        ▼                      data/filtered_bam_avg_quality.tsv
data/search_metadata.csv                   │
        │                                  │
        ▼                                  │
   run_covar.sh (runs the `covar` CLI)     │
        ▼                                  │
data/covar_out/*.covariants.tsv            │
        │  query_clinical_api.py           │
        │  (queries outbreak.info)         │
        ▼                                  │
data/covar_clinical_detections.tsv ◄────┐  │
        │                               │  │
        ▼                               │  │
snp_frequencies/snp_frequencies.ipynb   │  │
        ▼                               │  │
snp_frequencies/{cryptic,non_cryptic}_snp_freqs.csv
        │                               │  │
        ├──► fig3/snp-error-model/error-model.ipynb ──► iontorrent/illumina_error_matrix.npy
        │                                                        │
        └──► supp_fig3/depth-quality-model.ipynb ◄───────────────┘ (+ median_coverage_summary.tsv,
                        │                                             filtered_bam_avg_quality.tsv)
                        ▼
              data/scaling_factor.csv
                        │
                        ▼
      fig2/cryptics-linechart.ipynb, fig2/cryptics-heatmap.ipynb  (normalization)
```

`data/covar_clinical_detections.tsv` (one row per linked-mutation cluster per sample, with `num_clinical_detections` = global clinical prevalence from outbreak.info) is the file to start from if you just want to reproduce a figure without re-running alignment/coVar/API queries — it, along with the other small derived CSVs listed below, is already committed to the repo.

## Directory guide

### `data/`
Raw metadata and the aggregated tables everything else depends on.

- `scripts/` — the upstream pipeline, in run order: `download_search_data.sh` → `rename_reference.sh` → `parse_search_metadata.py` (writes `search_metadata.csv`) → `run_covar.sh` (writes `covar_out/*.covariants.tsv`) → `query_clinical_api.py` (writes `covar_clinical_detections.tsv`). `assign_bg_lineage.py` is an optional lineage-assignment step whose output isn't used by any committed figure (the heatmap notebook does its own lineage assignment inline instead). `run_ivar.sh` is a separate ivar-based variant-calling pass used only to build the sequencing-error background model (see `fig3/snp-error-model/`); `process_ngs.sh` is an unrelated/legacy alignment script for a different sample batch.

  **Note:** `query_clinical_api.py` currently queries the [outbreak.info](https://outbreak.info) mutation-prevalence API to get `num_clinical_detections` per cluster. This will need to be updated once the new Muninn SARS-CoV-2 clinical API goes live.
- `covar_out/` — per-sample coVar output (linked mutation clusters + depth/frequency). Only consumed in bulk by `query_clinical_api.py`.
- `covar_clinical_detections.tsv` — the aggregated, clinically-annotated table read directly by nearly every figure notebook (fig2, fig3, `snp_frequencies/`).
- `search_metadata.csv` — per-sample site/date metadata (parsed from BAM filenames). Used by `query_clinical_api.py`, `supp_fig2/plot_metadata.ipynb`, `supp_fig3/depth-quality-model.ipynb`, `fig3/snp-error-model/error-model.ipynb`.
- `sars2_metadata/` — reference genome/annotation (`NC_045512_Hu-1.*`, older `2019-nCoV.*`), primer scheme (`ARTICv5.3.2.bed`), Freyja `usher_barcodes.feather` and `lineages.yml` (used for lineage assignment in fig2's heatmap and fig3's SNP linechart script).
- `lineage-prevalence/`, `qPCR/` — per-WWTP (Encina, Point Loma, South Bay) Freyja lineage prevalence and qPCR viral-load time series, feeding `fig2/ww_plot.py` and `fig2/cryptics-linechart.ipynb`.
- `filtered_bam_avg_quality.tsv`, `scaling_factor.csv` — intermediate outputs of the depth/error-rate normalization model (see `supp_fig3/`), consumed by fig2's normalization steps.
- `excluded_samples.csv`, `search_bucket.txt`, `plot_config.yml` — supporting config for the pipeline/plots.

### `snp_frequencies/`
Builds the weekly SNP-substitution-type tables used by the sequencing-error model. `snp_frequencies.ipynb` reads `data/covar_clinical_detections.tsv` and bins single-nucleotide substitutions by week and type, writing `non_cryptic_snp_freqs.csv`. `cryptic_snp_freqs.csv` is the equivalent table restricted to cryptic clusters (`num_clinical_detections == 0`); as committed, the notebook itself only shows the non-cryptic code path, so `cryptic_snp_freqs.csv` should be treated as a fixed input generated by a filtered variant of that notebook rather than something you can currently regenerate by re-running the file as-is.

These two CSVs feed **two** downstream analyses:
1. `fig3/snp-error-model/error-model.ipynb` — combines `cryptic_snp_freqs.csv` with platform-specific background substitution-error matrices (built from `run_ivar.sh` output on intergenic regions) to error-correct the observed cryptic SNP spectrum (the Supp. Fig. 4 analysis), and to write `iontorrent_error_matrix.npy` / `illumina_error_matrix.npy`, which are then used by `fig3/scripts/evol_trajectory.py` to score/rank individual cryptic clusters for the Fig. 3C descent plots.
2. `supp_fig3/depth-quality-model.ipynb` — sums `cryptic_snp_freqs.csv` to weekly totals and regresses them against sequencing depth × error rate to derive `data/scaling_factor.csv`, the normalization factor used in Fig. 2.

`non_cryptic_snp_freqs.csv` is not consumed by any current figure script (it's referenced only in a commented-out line in `error-model.ipynb`), so it's effectively a background/comparison dataset kept for reference.

### `fig2/` — cryptic variant burden vs. lineage prevalence and viral load (Fig. 2)
- `ww_plot.py` — smooths per-WWTP Freyja lineage prevalence and qPCR viral load, writes `data/point_loma_prevalence_smoothed.csv` and the per-site deconvolution plots in `plots/`. (`ww_all_plot.py` is an earlier draft of the same logic, superseded by `ww_plot.py`.)
- `get_median_coverage.sh` — computes per-BAM median depth from `bam/`, writes `median_coverage_summary.tsv` (used by `supp_fig3/depth-quality-model.ipynb`).
- `cryptics-linechart.ipynb` (Fig. 2A/B) — bins cryptic-cluster counts by week from `data/covar_clinical_detections.tsv`, normalizes using `data/scaling_factor.csv`, and overlays smoothed viral load and lineage prevalence from `point_loma_prevalence_smoothed.csv`.
- `cryptics-heatmap.ipynb` (Fig. 2C) — assigns cryptic clusters to parent VOC lineages via `usher_barcodes.feather`/`lineages.yml`, normalizes by `scaling_factor.csv`, plots a month × lineage heatmap.
- `count_cryptics.py` — quick script to print summary counts of cryptic clusters (used to derive numbers quoted in text, not a figure panel).
- `plots/` — output PDFs for this figure.

### `fig3/` — persistence, diversification, and stepwise evolution of cryptic variants (Fig. 3)
- `distribution.ipynb` (Fig. 3A/B) — histograms of cryptic cluster detection counts and detection-span days, from `data/covar_clinical_detections.tsv` restricted to the spike gene.
- `scripts/evol_history.py` — builds parent/child (subset-of-mutations) descent DAGs across all recurring cryptic clusters genome-wide; the hand-picked outputs are copied into `selected_evoplots/`.
- `scripts/evol_trajectory.py` — the Fig. 3C generator: for specific seed clusters, finds descendant clusters, scores each SNP using the platform-appropriate error matrix from `snp-error-model/`, and renders the ranked descent DAGs and swarmplots into `descent_plots/`.
- `scripts/snp_frequencies_linechart.py` — standalone diagnostic plot of monthly SNP-type counts genome-wide; not part of a committed figure panel.
- `snp-error-model/error-model.ipynb` — builds the IonTorrent/Illumina background substitution-error matrices and the error-corrected SNP-frequency plots (Supp. Fig. 4); see `snp_frequencies/` above for its inputs. `score-cryptics.ipynb` is an earlier, incomplete/broken attempt at this scoring (superseded by the scoring logic in `evol_trajectory.py`) and `coverage_by_week.tsv` is an orphaned intermediate with no current producer/consumer — both kept for reference only.
- `diversification.ipynb` — counts distinct cryptic clusters emerging after each VOC's introduction; note its output path (`variant_introduction/`) doesn't exist in this checkout, so it needs that directory created before it can save.
- `descent_plots/`, `selected_evoplots/` — committed Fig. 3C/3D panel outputs.

### `supp_fig1/` — coVar vs. SAM_Refiner vs. Crykey benchmark (Supp. Fig. 1)
- `bygul/` — simulated read generation: `simulate_ngs.sh` uses the `bygul` simulator plus a JN.1.1 reference genome (`EPI_ISL_18967918.fasta`) to produce synthetic amplicon reads at controlled depths, which are then aligned/trimmed/variant-called as ground truth.
- `scripts/test_accuracy*.sh`, `benchmark_*.sh` — run coVar, SAM_Refiner, and Crykey (via the top-level `crykey_wastewater.py`) on the simulated data and benchmark runtime with `hyperfine`.
- `scripts/plot_accuracy.ipynb`, `plot_performance.ipynb` — read the `results/` outputs of the scripts above and produce the accuracy Venn diagram and runtime comparison plots.
- `environment.yml` — conda environment for the benchmark tools (does not include coVar itself, which is assumed pre-installed).

### `supp_fig2/` — sample availability over time (Supp. Fig. 2)
`plot_metadata.ipynb` bins `data/search_metadata.csv` by month/WWTP for the wastewater sample panel, and `sd_clinical_samples.tsv` (a static, externally-sourced table) for the clinical sample panel.

### `supp_fig3/` — depth/error-rate normalization model (Supp. Fig. 3)
`avg_quality.sh` computes per-BAM average base quality (writes `data/filtered_bam_avg_quality.tsv`). `depth-quality-model.ipynb` combines that with `fig2/median_coverage_summary.tsv`, `data/search_metadata.csv`, and `snp_frequencies/cryptic_snp_freqs.csv` to fit the depth×error-rate regression, producing `data/scaling_factor.csv` (used in Fig. 2) and the diagnostic plots in `plots/`.

## Notes on raw data

Raw BAMs (`bam/`, `data/filtered_bam/`, `data/original_bam/`) and ivar outputs (`data/ivar/`) are gitignored and not included in this repo; the `data/scripts/` pipeline documented above regenerates them from S3 given `data/search_bucket.txt`. All CSV/TSV intermediates needed to reproduce the figures from `data/covar_clinical_detections.tsv` onward are committed.
