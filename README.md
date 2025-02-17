# deseq2 project

A project with the code to rin the differential expression analysis

Run from the existing directory of the project:

``` bash
Rscript run.R --run <path/to/the/config.yaml>
```

## Dependencies

1.DESeq2 2.GGally 3.ggrepel 4.tidyverse 5.optparse 6.patchwork 7.yaml 8.data.table

## Sample data tables

### Metadata

| condition         | sample_id | experiment_condition | Day | plate              | well | ngs_id   |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| Target\_\_plate10 | SAMP10    | AP                   | 10  | plate100_0 ug.ml_1 | 9    | NGSXXXXX |

### targeting entity data(tar_entities.csv)

| id   | name | parent_id | gene  |
|------|------|-----------|-------|
| 1589 | 1589 | 1589      | SMAD4 |

### eff_entities.csv

| id   | entity_id | parent_id | name  |
|------|------|-----------|-------|
| 1589 | 1589 | 1589      | 1589 |

### gene_counts.tsv and gene_tpm.txv

These files have the count and tpm information

### cond_pairs.csv

File containing conditions and reference conditions to provide it to deseq wrapper.


