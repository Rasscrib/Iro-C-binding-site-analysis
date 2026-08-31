# Iro-C binding site analysis

R script for the *in silico* identification and visualization of Iro-C binding motifs
(ACAnnTGT) across the *upd1*/*upd2*/*upd3* locus, used in:

**Nakanishi T, Enomoto M, Igaki T.** *Iro-C/IRX creates anti-cancerized epithelial field against IL-6-dependent malignant tumorigenesis.*

`Iro-C_binding_site_plot.R` reproduces **Figure EV4G**.

## Requirements

R (tested with 4.4.3), and the packages `Gviz`, `BSgenome`,
`BSgenome.Dmelanogaster.UCSC.dm6`, `Biostrings`, `GenomicRanges`, `biomaRt`, and `dplyr`.

## Usage

Open `Iro-C_binding_site_plot.R` in R and run it. An internet connection is required
(the script retrieves gene models from Ensembl).

## Session information

```
R 4.4.3; Gviz 1.50.0, biomaRt 2.62.1, Biostrings 2.74.1, BSgenome 1.74.0,
BSgenome.Dmelanogaster.UCSC.dm6 1.4.1, GenomicRanges 1.58.0, dplyr 1.1.4
```
