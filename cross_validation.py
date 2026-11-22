"""
cross_validation.py
-------------------
K-fold cross-validation split for the Fair MFVI pipeline.

The function ``train_test_split_cv`` mirrors ``Code_2_fun_CV.m``.

For a given fold *n* it:
  1. Defines train / test index sets.
  2. Builds the fairness direction vector z for the *train* set and
     z_T for the *test* set following Eq. (23) of the paper:
         z = mean_features(non-sensitive, train) - mean_features(sensitive, train)
     Both z and z_T include a leading zero (for the intercept column) and
     are normalised to unit norm.

Public API
----------
train_test_split_cv(X0, y0, is_not_sensitive, M, n)
    -> (X, y, z, XT, yT, zT,
        ind_train, ind_test,
        mean_ns_train, mean_s_train,
        mean_ns_test,  mean_s_test)
"""

import numpy as np


def train_test_split_cv(
    X0: np.ndarray,
    y0: np.ndarray,
    is_not_sensitive: np.ndarray,
    M: int,
    n: int,
):
    """
    Create a single train/test split for fold *n* of a K-fold CV scheme.

    Parameters
    ----------
    X0 : ndarray, shape (N, D)
        Full design matrix (with intercept column).
    y0 : ndarray, shape (N,)
        Labels in {-1, +1}.
    is_not_sensitive : ndarray of bool, shape (N,)
        True for non-sensitive individuals, False for sensitive ones.
    M : int
        Approximate test-set size per fold  (= floor(N / num_folds)).
    n : int
        Fold index, 1-based  (n = 1, 2, ..., num_folds).

    Returns
    -------
    X : ndarray, shape (N_train, D)   -- train design matrix
    y : ndarray, shape (N_train,)     -- train labels
    z : ndarray, shape (D,)           -- fairness vector for train set
    XT: ndarray, shape (N_test,  D)   -- test  design matrix
    yT: ndarray, shape (N_test,)      -- test  labels
    zT: ndarray, shape (D,)           -- fairness vector for test  set
    ind_train : ndarray of int        -- indices of training samples
    ind_test  : ndarray of int        -- indices of test samples
    mean_ns_train : ndarray, shape (D-1,)
    mean_s_train  : ndarray, shape (D-1,)
    mean_ns_test  : ndarray, shape (D-1,)
    mean_s_test   : ndarray, shape (D-1,)
    """
    N = X0.shape[0]

    # ------------------------------------------------------------------
    # Build index sets using a boolean mask (avoids dtype issues with
    # list comprehensions over numpy integers)
    # Test  : [ (n-1)*M , n*M )  — 0-based, contiguous block
    # Train : everything else
    # ------------------------------------------------------------------
    test_start = int((n - 1) * M)
    test_end   = int(n * M)

    test_mask  = np.zeros(N, dtype=bool)
    test_mask[test_start:test_end] = True

    ind_test  = np.where( test_mask)[0]   # integer indices of test  samples
    ind_train = np.where(~test_mask)[0]   # integer indices of train samples

    X  = X0[ind_train, :]
    XT = X0[ind_test,  :]
    y  = y0[ind_train]
    yT = y0[ind_test]

    # ------------------------------------------------------------------
    # Fairness vectors  z  and  z_T
    # (Eq. 23)  z = mean_features(S=1) - mean_features(S=0),  ||z|| = 1
    # Features are all columns EXCEPT the intercept (column 0).
    # ------------------------------------------------------------------

    # -- Train --
    feats_train      = X0[ind_train, 1:]               # (N_train, D-1)
    sens_mask_train  = ~is_not_sensitive[ind_train]    # True -> sensitive

    feat_ns_train = feats_train[~sens_mask_train]      # non-sensitive rows
    feat_s_train  = feats_train[ sens_mask_train]      # sensitive rows

    mean_ns_train = feat_ns_train.mean(axis=0)         # (D-1,)
    mean_s_train  = feat_s_train.mean(axis=0)          # (D-1,)

    z_raw  = mean_ns_train - mean_s_train
    z_norm = np.linalg.norm(z_raw)
    z_dir  = z_raw / z_norm if z_norm > 0 else z_raw

    # Prepend zero for the intercept dimension
    z = np.concatenate([[0.0], z_dir])                 # (D,)

    # -- Test --
    feats_test      = X0[ind_test, 1:]
    sens_mask_test  = ~is_not_sensitive[ind_test]

    feat_ns_test = feats_test[~sens_mask_test]
    feat_s_test  = feats_test[ sens_mask_test]

    mean_ns_test = feat_ns_test.mean(axis=0)
    mean_s_test  = feat_s_test.mean(axis=0)

    zT_raw  = mean_ns_test - mean_s_test
    zT_norm = np.linalg.norm(zT_raw)
    zT_dir  = zT_raw / zT_norm if zT_norm > 0 else zT_raw

    zT = np.concatenate([[0.0], zT_dir])               # (D,)

    return (
        X, y, z,
        XT, yT, zT,
        ind_train, ind_test,
        mean_ns_train, mean_s_train,
        mean_ns_test,  mean_s_test,
    )