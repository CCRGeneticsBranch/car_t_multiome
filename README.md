Pipeline to process Illumina multimodal 10X Genomics single cell data, including:
- 5' GEX
- CITEseq (antibody capture)
- HTO multiplexing
- VDJ enrichment

Steps include:
1. Link fastq with bcl2fastq format names
2. Make cellranger multi config file
3. Run cellranger multi, including HTO demultiplexing
4. Perform alternative demultiplexing with MULTIseqDemux
5. Perform GEX and antibody normalization and filtering

To run the pipeline, clone the repo from github. In an interactive session, initialize the analysis directory with:
```
  bash [path to repo]/car_t_multiome/run_citeseq --runmode=init --workdir=[directory]
```
Edit the config/samples.tsv and config/config.yaml files as needed. Perform a dryrun with:
```
  bash [path to repo]/car_t_multiome/run_citeseq --runmode=dryrun --workdir=[directory]
```
If successful, launch the pipeline with:
```
  bash [path to repo]/car_t_multiome/run_citeseq --runmode=run --workdir=[directory]
```
