#!/usr/bin/env Rscript

# an example plotting script for use with the phyloseq R package.
# it is not intended to be a comprehensive analysis, but rather a starting point.
# for a tutorial on phyloseq, see: https://joey711.github.io/phyloseq/

# argument parsing
args = commandArgs(trailingOnly=TRUE)
biom_fp = args[1]
tax_fp = args[2]
metadata_fp = args[3]
output_fp = args[4]

# check for required packages
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package not found. please install via BiocManager.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package not found. please install.")
}
if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("biomformat package not found. please install via BiocManager.")
}

# load packages
library(phyloseq)
library(ggplot2)
library(biomformat)

# import qiime2 artifacts
biom_data <- import_biom(biom_fp)
tax_data <- read.table(tax_fp, sep="\t", header=TRUE, row.names=1)
metadata <- read.table(metadata_fp, sep="\t", header=TRUE, row.names=1)

# create phyloseq object
physeq <- phyloseq(otu_table(biom_data), 
                   tax_table(as.matrix(tax_data)), 
                   sample_data(metadata))

# create relative abundance plot
# this is a very basic example; for a real analysis, you would want to
# filter low-abundance taxa, agglomerate to a higher taxonomic level, etc.
p <- plot_bar(physeq, fill="Phylum") + 
    geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
    labs(x = "Sample", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

# save plot
ggsave(output_fp, plot=p, device="pdf", width=11, height=8.5)

print(paste("Plot saved to:", output_fp)) 