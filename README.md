# Iro-C binding site analysis

This repository contains the R script used for the *in silico* identification and
visualization of Iro-C binding motifs in the study:

**Nakanishi T, Enomoto M, Igaki T.**
*Iro-C/IRX creates anti-cancerized epithelial field against IL-6-dependent malignant tumorigenesis.*

The script reproduces the physical map of the *upd1*/*upd2*/*upd3* locus with Iro-C
consensus motifs shown in **Figure EV4G**.

## Contents

| File | Description |
|------|-------------|
| `Iro-C_binding_site_plot.R` | Retrieves the *upd1/upd2/upd3* locus sequence from the *D. melanogaster* reference genome (dm6), searches for the Iro-C minimal binding motif (ACAnnTGT) on the plus strand, retrieves gene models from Ensembl, and draws the locus map with Gviz. |

## What the script does

1. Obtains the genomic sequence spanning the *upd1*, *upd2* and *upd3* locus
   (chrX:18,230,000–18,313,000) from the *Drosophila melanogaster* reference genome
   (assembly dm6) via `BSgenome.Dmelanogaster.UCSC.dm6`.
2. Identifies occurrences of the Iro-C minimal binding motif **ACAnnTGT** (where n is any
   nucleotide) on the **plus strand** by regular-expression matching (`ACA[ACGT]{2}TGT`)
   using `Biostrings`.
3. Retrieves gene models within the region from Ensembl using `biomaRt`
   (dataset: `dmelanogaster_gene_ensembl`), keeping the longest transcript of each gene.
4. Visualizes the locus map, gene models, and motif positions using `Gviz`.

## Requirements

- **R** (tested with version 4.4.3)
- Bioconductor packages: `Gviz`, `BSgenome`, `BSgenome.Dmelanogaster.UCSC.dm6`,
  `Biostrings`, `GenomicRanges`, `biomaRt`
- CRAN package: `dplyr`

Install (once):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Gviz", "BSgenome.Dmelanogaster.UCSC.dm6",
                       "Biostrings", "GenomicRanges", "biomaRt"))
install.packages("dplyr")
```

## Usage

Open `Iro-C_binding_site_plot.R` in R (or RStudio) and run it. The script connects to
Ensembl to retrieve gene models, so an internet connection is required. If the Ensembl
connection times out, use a mirror by editing the `useEnsembl()` call, e.g.:

```r
mart <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl", mirror = "asia")
```

Running the script prints the number of identified motifs and draws the locus map
corresponding to Figure EV4G.

## Session information

The analysis was performed under the following environment:

```
R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64

other attached packages:
  dplyr 1.1.4
  biomaRt 2.62.1
  BSgenome.Dmelanogaster.UCSC.dm6 1.4.1
  BSgenome 1.74.0
  Biostrings 2.74.1
  Gviz 1.50.0
  GenomicRanges 1.58.0
  GenomeInfoDb 1.42.3
  IRanges 2.40.1
  S4Vectors 0.44.0
  BiocGenerics 0.52.0
```

## Citation

If you use this code, please cite the paper above.
