#!/usr/bin/env Rscript
# Wilcoxon rank-sum test on log2 TMM-normalised CPM, one test per gene.
# Genes are filtered with edgeR::filterByExpr, as for edgeR and limma-voom.
# logFC is the difference in mean log2-CPM (prior count 2) between snf2 and WT.
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(f[1])), "R", "common.R"))
})
suppressPackageStartupMessages(library(edgeR))

opt <- parse_de_args("wilcoxon")
inp <- read_de_inputs(opt)

runtime <- system.time({
  y <- DGEList(inp$counts, group = inp$coldata$condition)
  keep <- filterByExpr(y, design = inp$design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 2)

  is_ko <- inp$coldata$condition == "snf2"
  pvalue <- apply(lcpm, 1, function(x) {
    # Exact test where possible; R falls back to the normal approximation with ties.
    suppressWarnings(wilcox.test(x[is_ko], x[!is_ko])$p.value)
  })
  logFC <- rowMeans(lcpm[, is_ko, drop = FALSE]) - rowMeans(lcpm[, !is_ko, drop = FALSE])
})

res <- data.frame(gene = rownames(lcpm), logFC = logFC, pvalue = pvalue,
                  padj = p.adjust(pvalue, method = "BH"))
write_de_outputs(opt, inp, res, "Wilcoxon", c("edgeR", "stats"), runtime,
                 notes = "filterByExpr; TMM log2-CPM (prior.count = 2); stats::wilcox.test defaults")
