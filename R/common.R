# Shared input/output handling for the DE entrypoints.
#
# Every entrypoint receives raw counts and writes:
#   de_results.tsv  gene, logFC, pvalue, padj  (all input genes, input order;
#                   genes a method did not test have NA pvalue/padj)
#   de_info.json    method, package, versions, runtime, gene and sample counts
#
# logFC is log2(snf2 / WT). padj is Benjamini-Hochberg over the tested genes.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

parse_de_args <- function(method) {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--sub.counts", dest = "counts", type = "character"),
    make_option("--sub.coldata", dest = "coldata", type = "character"),
    make_option("--output_dir", type = "character"),
    make_option("--name", type = "character", default = method)
  )))
  for (a in c("counts", "coldata", "output_dir")) {
    if (is.null(opt[[a]])) stop("Missing required argument: ", a)
  }
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  opt
}

read_de_inputs <- function(opt) {
  tab <- read.delim(gzfile(opt$counts), check.names = FALSE)
  stopifnot(identical(colnames(tab)[1], "gene"))
  counts <- as.matrix(tab[, -1, drop = FALSE])
  rownames(counts) <- tab$gene
  storage.mode(counts) <- "integer"

  coldata <- read.delim(opt$coldata, colClasses = c(sample = "character", condition = "character"))
  coldata <- coldata[match(colnames(counts), coldata$sample), , drop = FALSE]
  if (anyNA(coldata$sample)) stop("Count columns and coldata samples do not match")
  if (!all(coldata$condition %in% c("WT", "snf2"))) stop("Unexpected condition labels")
  coldata$condition <- factor(coldata$condition, levels = c("WT", "snf2"))
  rownames(coldata) <- coldata$sample
  if (any(table(coldata$condition) < 1)) stop("Each condition needs at least one sample")

  list(counts = counts, coldata = coldata,
       design = model.matrix(~ condition, data = coldata))
}

# res: data.frame with columns gene, logFC, pvalue, padj (any subset of genes)
write_de_outputs <- function(opt, inputs, res, method, packages, runtime, notes = NULL) {
  out <- data.frame(gene = rownames(inputs$counts))
  res <- res[match(out$gene, res$gene), , drop = FALSE]
  out$logFC <- res$logFC
  out$pvalue <- res$pvalue
  out$padj <- res$padj
  write.table(out, file.path(opt$output_dir, "de_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  versions <- lapply(setNames(packages, packages), function(p) as.character(packageVersion(p)))
  info <- list(
    method = method,
    module = opt$name,
    packages = versions,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    runtime_seconds = round(unname(runtime[["elapsed"]]), 3),
    n_genes = nrow(out),
    n_tested = sum(!is.na(out$pvalue)),
    n_samples = as.list(table(inputs$coldata$condition)),
    notes = notes
  )
  write_json(info, file.path(opt$output_dir, "de_info.json"),
             auto_unbox = TRUE, pretty = TRUE, null = "null")
  message(sprintf("%s: %d/%d genes tested, %d with padj < 0.05", method,
                  info$n_tested, info$n_genes, sum(out$padj < 0.05, na.rm = TRUE)))
}
