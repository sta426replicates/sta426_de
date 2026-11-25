#!/usr/bin/env Rscript
# limma-voom: filterByExpr, TMM normalisation, voom, empirical Bayes moderated t-test.
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(f[1])), "R", "common.R"))
})
suppressPackageStartupMessages({ library(edgeR); library(limma) })

opt <- parse_de_args("limma_voom")
inp <- read_de_inputs(opt)

runtime <- system.time({
  y <- DGEList(inp$counts, group = inp$coldata$condition)
  keep <- filterByExpr(y, design = inp$design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, inp$design)
  fit <- eBayes(lmFit(v, inp$design))
  tt <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
})

res <- data.frame(gene = rownames(tt), logFC = tt$logFC, pvalue = tt$P.Value, padj = tt$adj.P.Val)
write_de_outputs(opt, inp, res, "limma-voom", c("limma", "edgeR"), runtime,
                 notes = "filterByExpr; TMM; voom; eBayes")
