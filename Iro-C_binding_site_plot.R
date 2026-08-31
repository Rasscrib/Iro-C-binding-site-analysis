############################################################
## Identification and visualization of Iro-C binding motifs
## (ACAnnTGT) across the upd1/upd2/upd3 locus
##
## Nakanishi et al., "Iro-C/IRX creates anti-cancerized
## epithelial field against IL-6-dependent malignant
## tumorigenesis"
##
## Genome : D. melanogaster dm6 (BSgenome.Dmelanogaster.UCSC.dm6)
## Region : chrX:18,230,000-18,313,000 (upd1/upd2/upd3 locus)
## Motif  : ACAnnTGT (Iro-C minimal binding site), plus strand
## Output : locus map (axis, gene models, Iro-C motifs)
##
## This single script reproduces Figure EV4G.
############################################################

## ---- Installation (run once; uncomment if packages are not installed) ----
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("Gviz", "BSgenome.Dmelanogaster.UCSC.dm6",
#                        "Biostrings", "GenomicRanges", "biomaRt"))
# install.packages("dplyr")

## ---- Load packages (all required libraries are attached here) ----
library(BSgenome)                        # provides getSeq()
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(Biostrings)
library(GenomicRanges)
library(biomaRt)                         # provides useEnsembl() / getBM()
library(dplyr)                           # provides %>%, group_by(), filter()
library(Gviz)

## ---- Target region ----
chr  <- "chrX"
from <- 18230000L
to   <- 18313000L

## ---- Iro-C motif (ACAnnTGT, plus strand only) ----
s_obj <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, chr, start = from, end = to)
s     <- if (is(s_obj, "DNAStringSet")) s_obj[[1]] else s_obj
ss    <- as.character(s)

pattern_regex <- "ACA[ACGT]{2}TGT"       # ACAnnTGT (8 bp)

m_f      <- gregexpr(pattern_regex, ss, perl = TRUE, ignore.case = TRUE)
starts_f <- as.integer(m_f[[1]]); starts_f <- starts_f[starts_f > 0]
ends_f   <- if (length(starts_f) > 0) starts_f + 8L - 1L else integer(0)

# convert local coordinates to genomic coordinates
start_genome_f <- if (length(starts_f) > 0) from + starts_f - 1L else integer(0)
end_genome_f   <- if (length(ends_f)   > 0) from + ends_f   - 1L else integer(0)

cat("Number of Iro-C motif hits (plus strand):", length(starts_f), "\n")

motifTrack <- AnnotationTrack(
  start      = start_genome_f, end = end_genome_f,
  chromosome = chr, genome = "dm6",
  name       = "Iro-C motif (+)",
  shape      = "box", fill = "#e41a1c", col = "#e41a1c",
  size       = 0.20, stacking = "dense"
)

## ---- Gene models (longest transcript per gene, grey, italic labels) ----
mart <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")
# If the connection fails or times out, specify a mirror, e.g.:
# mart <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl", mirror = "asia")

# Retrieve all exons and transcript lengths within the region
exon_df <- getBM(
  attributes = c("external_gene_name", "ensembl_transcript_id", "chromosome_name",
                 "strand", "exon_chrom_start", "exon_chrom_end", "rank",
                 "transcript_length"),
  filters    = c("chromosome_name", "start", "end"),
  values     = list("X", from, to),
  mart       = mart
)
exon_df$chromosome_name <- paste0("chr", exon_df$chromosome_name)

# Keep only the longest transcript for each gene
# (dplyr::filter is namespaced to avoid masking by Biostrings/S4Vectors::filter)
exon_df_longest <- exon_df %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::filter(transcript_length == max(transcript_length, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Build a GRanges object of exons
exon_gr <- GRanges(
  seqnames   = exon_df_longest$chromosome_name,
  ranges     = IRanges(exon_df_longest$exon_chrom_start,
                       exon_df_longest$exon_chrom_end),
  strand     = ifelse(exon_df_longest$strand == 1, "+", "-"),
  transcript = exon_df_longest$ensembl_transcript_id,
  symbol     = exon_df_longest$external_gene_name,
  exon       = exon_df_longest$rank
)

geneTrack <- GeneRegionTrack(
  range                = exon_gr,
  genome               = "dm6", chromosome = chr,
  transcriptAnnotation = "symbol",    # label with gene symbol
  groupAnnotation      = "symbol",    # group by gene
  collapseTranscripts  = "longest",   # collapse each gene to its longest transcript
  fill = "grey60", col = "grey30",    # all grey
  fontcolor.group = "black",
  name = "Genes"
)

# gene labels in italic, placed to the left
displayPars(geneTrack) <- list(
  fontface.group = "italic",
  just.group     = "left"
)

## ---- Genome axis ----
axisTrack <- GenomeAxisTrack()

## ---- Plot (reproduces Figure S4G) ----
plotTracks(
  list(axisTrack, geneTrack, motifTrack), 
  from = from, to = to, chromosome = chr,
  betweenTrackMargin = 0,
  sizes    = c(1, 0.3, 0.3),   # relative track heights (adjust as needed)
  cex.title = 0.9, cex.axis = 0.8, cex.group = 0.85
)
