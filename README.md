# Short Reads vs Long Reads: settling it with data (E. coli K-12 MG1655)

Course project 18, *Introduction to Bioinformatics*, Astana IT University. Author: Otegen Symbat (individual project).

**Question.** For one well-characterised organism with a finished reference, what do short-read, long-read and hybrid
strategies each deliver in contiguity, base accuracy, repeat resolution and cost, and when is the long-read premium justified?

## Data (all real, public)
| Role | Accession | Notes |
|---|---|---|
| Reference | NC_000913.3 (RefSeq GCF_000005845.2) | E. coli K-12 MG1655, 4,641,652 bp |
| Short reads | SRR1030394 | Illumina MiSeq 2x251, ~285x |
| Long reads | ERR14686234 | Nanopore MinION, ~36x, **ciprofloxacin-resistant mutant C57** of MG1655 |

Downloads are scripted (`scripts/fetch_ena.py`, md5 verified, retrieval dates in `data/PROVENANCE.tsv`; the Colab run also saved it as `results/PROVENANCE.tsv` on Google Drive, copy it into this repository). No data files are committed.

## Reproduce (in Google Colab, how the results were produced)
1. Open `notebooks/01_pipeline_part1.ipynb` in Colab. Run cell 1 (it restarts the runtime), then run all remaining cells.
2. Open `notebooks/02_analysis_part2.ipynb` and run all cells. Results land in `results/`, large files stay on Google Drive.
3. `python scripts/make_figures.py` regenerates `figures/` from `results/*.csv`.

## Reproduce (Snakemake, **not yet run end to end**)
1. `mamba env create -f environment.yml && mamba activate p18`
2. `snakemake --cores 2 -n` (dry run), then `snakemake --cores 2`
(One command: `snakemake --cores 2` regenerates `results/quast/transposed_report.tsv`, `results/indel_context.csv` and the figures.)
Expect about 3 hours on 2 cores.

## Layout
- `notebooks/` executed notebooks (outputs included), `scripts/` standalone scripts, `Snakefile` + `config.yaml` workflow
- `results/` small result tables and QC summaries, `figures/` report figures, `tests/` tiny test data (see below)

## Status
| Item | State |
|---|---|
| Data retrieval, strain check, QC | done |
| Assemblies: short / long / hybrid (SPAdes hybrid, Flye + Polypolish) | done (Unicycler failed twice, see report) |
| QUAST evaluation, uncovered regions, homopolymer error analysis | done |
| Cost-normalised comparison | **TODO**: needs sourced price per Gb |
| BLAST of the extra ~85 kb contig | **TODO** |
| Snakemake workflow executed end to end | **TODO** |
| Small test dataset (`scripts/make_test_data.sh`) | **TODO**: script written, untested |
| BUSCO/annotation completeness, structural variants | not evaluated (stated as limitation) |
