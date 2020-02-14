"""
"
" Compute feature-phenotype spearman correlation matrix
"
" Spearman correlation for each feature (e.g. gene expression profile) and each 
" phenotype (e.g. drug response profile). The maximum correlation observed for each
" feature across all phenotypes is then returned.
"
"""
from numba import njit
import numpy as np

"""
Fast Spearman correlation calculation

https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes
"""
@njit
def mean1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].mean()
  return b

@njit
def std1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].std()
  return b

@njit
def spearman_cor(a, b):
  """
  Spearman correlation

  Expects value rankings from unrolled 2d arrays with the same numbers of columns
  """
  n, k = a.shape
  m, k = b.shape

  mu_a = mean1(a)
  mu_b = mean1(b)
  sig_a = std1(a)
  sig_b = std1(b)

  out = np.empty((n, m))

  for i in range(n):
    for j in range(m):
      out[i, j] = (a[i] - mu_a[i]) @ (b[j] - mu_b[j]) / k / sig_a[i] / sig_b[j]

  return out
