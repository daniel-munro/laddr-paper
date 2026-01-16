# laddr-paper

Code and outputs for the latent data-driven RNA phenotyping (LaDDR) paper. This repo holds the R scripts used to assemble intermediate datasets, compute paper statistics, and generate the figures and tables.

## Repository map

- `scripts/`: R scripts for assembling results, computing stats for the manuscript, and generating each figure/table.
- `data/`: inputs and intermediate data (QTLs, TWAS hits, phenotype groups, reference annotations, QC, etc.).
- `figures/`: final figure PNGs and design sources (`.afdesign`). R scripts in `scripts/` generate either individual figure panels or entire figures, and in the former case, the panel images are linked within an Affinity Designer file that lays out all assembled figures, which are exported as PNGs.
- `tables/`: supplementary tables (TSV and XLSX). R scripts in `scripts/` generate raw TSV tables, which I then copy into an Excel spreadsheet alongside a sheet with the title, column descriptions, etc.
- `laddr-paper.Rproj`: RStudio project file.

## Scripts at a glance

- `scripts/assemble_results.R`: builds processed datasets from raw inputs (genes, QTLs, TWAS hits, held-out analyses, pruning, sequencing sims, etc.).
- `scripts/stats.R`: computes the numeric values used in the paper text from processed data.
- `scripts/figure*.R`: generate main-text figures (e.g., `figure1a.R`, `figure2.R`, `figure3abc.R`, `figure4abcd.R`, `figure5a.R`).
- `scripts/figureS*.R`: generate supplementary figures.
- `scripts/tableS*.R`: generate supplementary tables.
- `scripts/latent_explicit_corrs.R`: analysis helper for latent vs explicit phenotype correlations.

## Data notes

- `data/download.sh` documents how some inputs (notably TWAS hits) were pulled via `rsync` from an internal server.
- Processed outputs referenced by `scripts/stats.R` and the figure/table scripts live under `data/processed/`.

## Outputs

- `figures/figure*.png` and `figures/figureS*.png` are the rendered figure outputs used in the manuscript.
- `tables/tableS*.tsv` and `tables/Supplementary_Data_*.xlsx` are the supplementary data files.
