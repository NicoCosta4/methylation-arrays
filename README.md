# methylation-arrays

Reproducible analysis of Illumina Infinium DNA methylation arrays (EPIC / EPICv2) with R and Bioconductor, packaged in a Docker container.

> **Status:** work in progress — the pipeline is being built step by step.

## Example dataset

The pipeline is demonstrated on [GSE245203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245203): EPIC (v1) methylation profiles of **3 primary tumors and 5 metastases** from a spontaneous metastasis mouse model of Ewing sarcoma (human A673 cells), published in:

> Chicon-Bosch M, *et al.* Multi-omics profiling reveals key factors involved in Ewing sarcoma metastasis. *Molecular Oncology* (2025). [doi:10.1002/1878-0261.13788](https://doi.org/10.1002/1878-0261.13788)

Metastases come from the same mouse as their primary tumor, so the samples are not independent. The differential methylation model accounts for this by treating the mouse as a blocking factor.

| Mouse | Primary tumor | Metastases |
|---|---|---|
| M6 | 1 | 0 |
| B3 | 1 | 2 |
| B5 | 1 | 3 |

## Pipeline

1. **Download** raw IDAT files from GEO
2. **Quality control** — detection p-values per sample and per probe
3. **Normalization** — stratified quantile normalization (`minfi::preprocessQuantile`)
4. **Probe filtering** — failed probes, sex chromosomes, SNP-affected and cross-reactive probes
5. **Differential methylation** — `limma` on M-values, mouse as blocking factor, FDR + Δβ thresholds
6. **Figures** — PCA, density plots, volcano plot and heatmap

## Usage

*Instructions will be added as the pipeline is completed.*

## License

[MIT](LICENSE)
