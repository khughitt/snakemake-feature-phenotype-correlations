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

_Example feature data_:

|symbol   |   1321N1|     143B|    22RV1| 23132-87|
|:--------|--------:|--------:|--------:|--------:|
|A1BG     | 5.542727| 5.245888| 4.753565| 4.069646|
|A1BG-AS1 | 4.658156| 4.560726| 4.766875| 4.471635|
|A1CF     | 3.921564| 3.776961| 7.445068| 6.272019|
|A2M      | 4.236747| 4.504467| 4.824879| 3.731857|
|A2M-AS1  | 3.535469| 4.376824| 5.714456| 4.850112|

_Example phenotype data_:

|drug      |    1321N1|     22RV1|  42-MG-BA|      5637|
|:---------|---------:|---------:|---------:|---------:|
|17-AAG    | 0.4177000| 0.3724600| 0.5990000| 0.4828250|
|AEW541    | 0.0873750| 0.2205000| 0.1144375| 0.1243550|
|AZD0530   | 0.2557625| 0.0262375| 0.1621000| 0.1627625|
|AZD6244   | 0.1107000| 0.3551250| 0.1045875| 0.0647625|
|Erlotinib | 0.0301375| 0.0256750| 0.0472750| 0.2240250|

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
