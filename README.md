# sta425_de

Omnibenchmark **de** module for the STA425 replicate benchmark: one repository,
one entrypoint per method.

| Entrypoint | Method | Filtering | Normalisation | Test |
|---|---|---|---|---|
| `deseq2` | DESeq2 | independent filtering (`alpha = 0.05`), Cook's outliers | median-of-ratios | Wald |
| `edger` | edgeR | `filterByExpr` | TMM | quasi-likelihood F-test |
| `limma_voom` | limma-voom | `filterByExpr` | TMM | moderated t (voom + eBayes) |
| `wilcoxon` | Wilcoxon rank-sum | `filterByExpr` | TMM log2-CPM | `wilcox.test` per gene |

All methods receive raw counts. The contrast is snf2 vs WT, so positive `logFC`
means higher in snf2.

## Interface

```
<entrypoint>.R --sub.counts counts.tsv.gz --sub.coldata coldata.tsv \
               --output_dir DIR --name ID
```

| Output | Format |
|---|---|
| `de_results.tsv` | `gene`, `logFC`, `pvalue`, `padj` — every input gene in input order; genes a method did not test have `NA` p-values |
| `de_info.json` | method, package versions, R version, runtime, number of genes tested, samples per condition |

`padj` is Benjamini–Hochberg over the genes each method tested. For DESeq2,
`logFC` is the unshrunken maximum-likelihood estimate.

## Environment

R with `DESeq2`, `edgeR`, `limma`, `optparse` and `jsonlite`
(`envs/r_de.yaml` in the benchmark repository).
