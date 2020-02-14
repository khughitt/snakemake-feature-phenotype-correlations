Snakemake Feature-Phenotype Correlation Pipeline
================================================

Overview
--------

A simple [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based pipeline for
generating correlation matrices within and across paired datasets (i.e. those with
measurements for a shared set of variables or samples).

Given a set of input feature and phenotype datasets, the pipeline constructs three types
of correlation matrices:

- **Feature correlation matrix** Row-wise feature correlations (e.g. gene-gene correlations)
- **Phenotype correlation matrix** Row-wise phenotype correlations (e.g. gene-gene correlations)
- **Feature-Phenotype correlation matrix** Feature-phenotype cross-dataset correlations for a set of shared columns (samples).

At present, the pipeline supports generating Pearson and Spearman corelation matrices.

Input Data
----------

The pipeline expects one or more "feature" datasets (e.g. mRNA expression measured at
the gene level), as well as one or more "phenotype" datasets (e.g. drug sensitivity
AUC scores), for a shared set of samples.

In each case, the datasets should be oriented with features/phenotypes along the rows,
and sample ids across columns. A header row should be included with the sample ids (it's
okay if they do not match exactly, as long as their is some subset of shared ids
present). The first column of each dataset should correspond to row identifiers. All
input data should be stored using the [feather format](https://github.com/wesm/feather).

Ex.:

```r
head(feats)

   symbol   1321N1     143B
     A1BG 5.542727 5.245888
 A1BG-AS1 4.658156 4.560726
     A1CF 3.921564 3.776961
      ...


head(phenos)

   drug    1321N1     22RV1
 17-AAG 0.4177000 0.3724600
 AEW541 0.0873750 0.2205000
AZD0530 0.2557625 0.0262375
    ...
```

Example Usage
-------------

To run the pipeline, copy and modify the example config file located in the `config/`
folder to point to the locations of your input datasets.

Next, create a conda environment and install the necessary required dependencies using
the following commands:


```sh
conda create -n fcor --file requirements.txt
```

You can then activate and run the pipeline with:

```sh
conda activate fcor
snakemake --configfile config.yml -j4
```
