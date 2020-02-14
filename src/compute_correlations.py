#!/bin/env python
"""
"
" Computes correlation matrix for columns of a single dataset
" KH Feb 2020
"
"""
import numpy as np
import pandas as pd

# load data
dat = pd.read_feather(snakemake.input[0])

# keep track of identifier field name
id_field = dat.columns[0]

# separate identifiers from data
row_ids = dat[id_field]
dat = dat.drop(id_field, axis=1)

# compute correlation matrix using specified method
if snakemake.wildcards["method"] == "pearson":
    # pearson correlation
    cor_mat = np.corrcoef(dat)
elif snakemake.wildcards["method"] == "spearman":
    # spearman correlation
    from scipy.stats import spearmanr

    cor_mat = spearmanr(dat.T).correlation

# convert result to a dataframe and update column and row names
cor_mat = pd.DataFrame(cor_mat)
cor_mat.columns = row_ids
cor_mat.index = row_ids

# move index to a column to allow saving to feather
cor_mat = cor_mat.reset_index().rename(columns={"index": id_field})

cor_mat.to_feather(snakemake.output[0])
