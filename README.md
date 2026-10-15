# Single-Cell Multiome Pipeline

Snakemake pipeline for multimodal Illumina 10X Genomics single-cell data: 5' gene expression (GEX), CITE-seq antibody capture, hashtag (HTO) sample multiplexing, and TCR (VDJ) sequencing (optional). CITE-seq measures cell-surface protein abundance using oligonucleotide-tagged antibodies. Cell hashing uses tagged antibodies against ubiquitous cell surface proteins to multiplex samples.

The pipeline runs Cell Ranger `multi` per library, demultiplexes hashtags, normalizes RNA and antibody data, applies QC filters, and saves Seurat objects per library and per sample. It supports single-cell multiome experiments across a range of research applications.

## References

The following references were used to build the pipeline: 

+ https://satijalab.org/seurat/articles/hashing_vignette
+ https://satijalab.org/seurat/articles/multimodal_vignette
+ https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
+ Code developed by Meijie Tian and Shaoli Das
+ [CCBR SINCLAIR](https://ccbr.github.io/SINCLAIR/)

## Pipeline steps

| Step | Rule | Tool | Description |
|---|---|---|---|
| 1 | `link_fastq` | | Link FASTQs into `results/raw/` with bcl2fastq format names |
| 2 | `make_multi_csv` | Python | Write the Cell Ranger `multi` config file per library |
| 3 | `cellranger` | Cell Ranger `multi` | Alignment, counting, VDJ assembly, and hashtag assignment per library |
| 4 | `copy_summary` | | Copy each sample's `web_summary.html` to `results/cellranger/` |
| 5 | `preprocess_sample` | R (Seurat, dsb) | From the Cell Ranger filtered sample objects, performs RNA normalization, DSB antibody normalization using empty-droplet background, and cell filtering |
| 6 | `preprocess_library` | R (Seurat, dsb) | From the Cell Ranger raw library object, demultiplexes hashtags with `MULTIseqDemux`, then performs normalization and filtering |

DSB normalization removes ambient noise using empty droplets and uses non-specific isotype controls to remove inter-cell technical variation. 

## Requirements

- Runs on NIH [Biowulf](https://hpc.nih.gov/) (SLURM and environment modules).
- Must have access to dependencies listed in [config/config.yaml](config/config.yaml).
- R packages used by the reports: Seurat, dsb, dplyr, ggplot2, patchwork, reshape2, ggridges, and rmarkdown.

## Installation

```bash
git clone https://github.com/CCRGeneticsBranch/car_t_multiome.git
```

## Quick start

The pipeline is designed to be run in a separate directory (/path/to/my_run).

```bash
# 1. Create a work folder with template config files
bash car_t_multiome/run_car_t_multiome --runmode=init --workdir=/path/to/my_run [--samplesheet=/path/to/samples.tsv]

# 2. Edit /path/to/my_run/config/config.yaml and config/samples.tsv (see Inputs)

# 3. Dry run (writes /path/to/my_run/dryrun.log)
bash car_t_multiome/run_car_t_multiome --runmode=dryrun --workdir=/path/to/my_run

# 4. Submit to SLURM
bash car_t_multiome/run_car_t_multiome --runmode=run --workdir=/path/to/my_run

# Unlock the work folder after an interrupted run
bash car_t_multiome/run_car_t_multiome --runmode=unlock --workdir=/path/to/my_run
```

`-h` shows all options; `-v` prints the pipeline and launcher versions. Per-rule SLURM resources are in `WORKDIR/config/cluster.yaml` (copied from [resources/cluster.yaml](resources/cluster.yaml) at `init`).

## Inputs

### Sample sheet (`config/samples.tsv`)

Tab-separated, one row per sample (hashtag) in each library. If a library was sequenced on several flowcells, each instance should be listed with the corresponding flowcell ID. 

| Column | Description | Example |
|---|---|---|
| `Library` | Library name; part of the FASTQ folder name | `library1` |
| `FCID` | Flowcell ID; part of the FASTQ folder name | `FLOWCELL1` |
| `Sample` | Sample name within the library | `sample1` |
| `Hashtag` | Hashtag ID for the sample; must match an ID in the feature reference CSV | `HTO1` |

### Raw FASTQ layout

FASTQs must be in `rawdir`, in one folder per library, mode, and flowcell:

```
rawdir/Sample_<Library>_<mode>_<FCID>/Sample_<Library>_<mode>_<FCID>_R1.fastq.gz
rawdir/Sample_<Library>_<mode>_<FCID>/Sample_<Library>_<mode>_<FCID>_R2.fastq.gz
```

`<mode>` is `5GEX` (gene expression), `BC` (antibody and hashtag capture), or `TCR` (VDJ sequencing).

### Configuration (`config/config.yaml`)

| Key | Default | Description (source) | Edit? |
|---|---|---|---|
| `rawdir` | `/data/khanlab/projects/DATA` | Folder containing `Sample_*` FASTQ folders | As needed |
| `gex_ref` | `/path/to/custom/reference/` | Cell Ranger GEX reference (standard or custom; see [10x reference downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads)) | Yes |
| `feature_ref` | `/path/to/feature_ref.csv` | Cell Ranger [feature reference CSV](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv) for antibody and hashtag barcodes | Yes |
| `vdj_ref` | `/fdb/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/` | Cell Ranger VDJ reference | As needed |
| `has_tcr` | `true` | Include TCR libraries | Set `false` for GEX + CITE-seq only |
| `force` / `force_cells` | `false` / `20000` | Force Cell Ranger to call a specified number of cells | As needed |
| `nFeature_RNA_min` / `nFeature_RNA_max` | `500` / `6000` | Retained range of detected genes per cell | As needed |
| `percent_mt_max` | `10` | Maximum mitochondrial read percentage per cell | As needed |
| `workdir`, `manifest`, `tools`, `cluster` | `WORKDIR/...`, `PIPELINE_HOME/...` | Work folder and launcher/repository paths (filled by `init`) | No |

## Outputs

All paths are relative to the work folder.

| Path | Description |
|---|---|
| `results/raw/` | Renamed links to the input FASTQs |
| `results/cellranger/multi_<library>.csv` | Cell Ranger `multi` config |
| `results/cellranger/<library>/` | Cell Ranger `multi` output (per-sample gene x cell matrices, VDJ results) |
| `results/cellranger/<library>.<sample>.web_summary.html` | Cell Ranger summary per sample |
| `results/preprocessing/by_library/<library>.rds` / `.html` | Demultiplexed, normalized, and filtered Seurat object per library, with report |
| `results/preprocessing/by_sample/<library>.<sample>.rds` / `.html` | Normalized and filtered Seurat object per sample, with report |
| `results/pipeline_version.txt` | Pipeline version used for the run |
| `runinfo.yaml`, `snakemake.log`, `logs/` | Run information, Snakemake log and per-rule SLURM logs |

## Tools and references

Modules are listed in [resources/tools.yaml](resources/tools.yaml).

| Tool | Module / source | Used in | Reference |
|---|---|---|---|
| Snakemake | `snakemake/7.19.1` | workflow | Mölder et al. 2021, *F1000Research* [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| Cell Ranger | `cellranger/9.0.1` | `cellranger` | Zheng et al. 2017, *Nat Commun* [doi:10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049); [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger) |
| Seurat | `R` | `preprocess_*` | Hao et al. 2024, *Nat Biotechnol* [doi:10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y) |
| MULTIseqDemux (Seurat) | `R` | `preprocess_library` | McGinnis et al. 2019, *Nat Methods* [doi:10.1038/s41592-019-0433-8](https://doi.org/10.1038/s41592-019-0433-8) |
| dsb | `R` | `preprocess_*` | Mulè et al. 2022, *Nat Commun* [doi:10.1038/s41467-022-29356-8](https://doi.org/10.1038/s41467-022-29356-8) |
| CITE-seq / cell hashing | — | — | Stoeckius et al. 2017, *Nat Methods* [doi:10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380); Stoeckius et al. 2018, *Genome Biol* [doi:10.1186/s13059-018-1603-1](https://doi.org/10.1186/s13059-018-1603-1) |

## Versioning

The pipeline version is in [VERSION](VERSION) and changes are listed in [CHANGELOG.md](CHANGELOG.md). Each run records the pipeline version, git commit and launcher version in `runinfo.yaml` and `results/pipeline_version.txt`. The launch script is a copy of [biowulf-snakemake-launcher](https://github.com/CCRGeneticsBranch/biowulf-snakemake-launcher) (`bin/launcher.sh`).

## Contact

Erica Pehrsson, Genetics Branch, Center for Cancer Research, National Cancer Institute. Please report problems through GitHub issues.
