TODO

When analyzing coordinated interactions for subset of variants which can be tissue-specific variants or other flavors
of subsets you can relay on the given scripts with some additional flags as will be demonstrated here.
Last step which is the coordination test itself is slightly different, and therefore we used a separate `R` script for
that matter.

1. Summarize over all tissues
```shell script
summarize_tissue_specific_PRS.py --phenotype $f --pvalue $p --AllPCs
```
1. Last step uses:
```shell script
prs_summary_tissue_specific.R
```