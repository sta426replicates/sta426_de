#!/usr/bin/env Rscript
# edgeR: filterByExpr, TMM normalisation, quasi-likelihood F-test.
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(f[1])), "R", "common.R"))
})
suppressPackageStartupMessages(library(edgeR))

opt <- parse_de_args("edger")
inp <- read_de_inputs(opt)

runtime <- system.time({
  y <- DGEList(inp$counts, group = inp$coldata$condition)
  keep <- filterByExpr(y, design = inp$design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, inp$design)
  fit <- glmQLFit(y, inp$design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")$table
})

res <- data.frame(gene = rownames(tt), logFC = tt$logFC, pvalue = tt$PValue, padj = tt$FDR)
write_de_outputs(opt, inp, res, "edgeR", c("edgeR", "limma"), runtime,
                 notes = "filterByExpr; TMM; glmQLFit + glmQLFTest")
