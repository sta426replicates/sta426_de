#!/usr/bin/env Rscript
# DESeq2: Wald test, median-of-ratios size factors, default independent filtering
# and Cook's distance outlier handling. logFC is the unshrunken MLE.
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(f[1])), "R", "common.R"))
})
suppressPackageStartupMessages(library(DESeq2))

opt <- parse_de_args("deseq2")
inp <- read_de_inputs(opt)

runtime <- system.time({
  dds <- DESeqDataSetFromMatrix(inp$counts, inp$coldata, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, contrast = c("condition", "snf2", "WT"), alpha = 0.05)
})

res <- data.frame(gene = rownames(r), logFC = r$log2FoldChange,
                  pvalue = r$pvalue, padj = r$padj)
write_de_outputs(opt, inp, res, "DESeq2", c("DESeq2"), runtime,
                 notes = "Wald test; independent filtering with alpha = 0.05; Cook's cutoff default")
