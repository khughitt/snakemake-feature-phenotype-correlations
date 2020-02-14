#!/bin/env python
"""
"
" Computes correlation matrix for columns across two datasets with observations for
# the same samples/columns.
#
" KH Feb 2020
"
"""
import numpy as np
import pandas as pd

# load processed feature / phenotype data
feat_dat = pd.read_feather(snakemake.input['features'])
phen_dat = pd.read_feather(snakemake.input['phenotypes'])

# separate data and row ids
feat_id_field = feat_dat.columns[0]
phen_id_field = phen_dat.columns[0]

# separate identifiers from data
feat_ids = feat_dat[feat_id_field]
phen_ids = phen_dat[phen_id_field]

feat_dat = feat_dat.drop(feat_id_field, axis=1)
phen_dat = phen_dat.drop(phen_id_field, axis=1)

# get columns/samples present in both the feature and phenotype data;
shared_ids = feat_dat.columns
shared_ids = sorted(list(set(phen_dat.columns).intersection(shared_ids)))

# drop any columns that aren't shared
feat_dat  = feat_dat.loc[:, shared_ids]
phen_dat = phen_dat.loc[:, shared_ids]

# generate numpy versions of the datasets and transpose
feat_mat = feat_dat.to_numpy().T
phen_mat = phen_dat.to_numpy().T

# compute feat-phenotype correlation matrix
if snakemake.wildcards["method"] == "pearson":
    # pearson correlation
    from methods.pearson_correlation import pearson_cor

    cor_mat = pearson_cor(feat_mat, phen_mat).T

elif snakemake.wildcards["method"] == "spearman":
    # spearman correlation
    from methods.spearman_correlation import spearman_cor

    cor_mat = spearman_cor(pd.DataFrame(feat_mat.T).rank(1).values, 
                           pd.DataFrame(phen_mat.T).rank(1).values)

# convert back to a dataframe and add row/column identifiers
cor_mat = pd.DataFrame(cor_mat, columns = phen_ids)
cor_mat.insert(0, feat_id_field, feat_ids)

# store correlation matrix
cor_mat.to_feather(snakemake.output[0])
