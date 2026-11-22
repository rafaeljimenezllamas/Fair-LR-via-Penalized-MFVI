"""
data_loader.py
--------------
Dataset loading and preprocessing for the Fair MFVI experiments.

Two datasets are supported:

  - German Credit  (``load_german``)
  - Student Performance / Portuguese (``load_student``)

Both functions return a triple  (X0, y0, is_not_sensitive)  with
consistent conventions used throughout the pipeline:

  X0              : ndarray, shape (N, D)
                    Design matrix with a leading column of ones (intercept).
  y0              : ndarray, shape (N,)
                    Labels in {-1, +1}.
  is_not_sensitive: ndarray of bool, shape (N,)
                    True  -> individual belongs to the NON-sensitive group.
                    False -> individual belongs to the sensitive group.
"""

import numpy as np
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# German Credit
# ---------------------------------------------------------------------------

def load_german(filepath: str | Path = "german_credit.csv"):
    """
    Load and preprocess the German Credit dataset.

    Sensitive variable: Age (threshold > 25 years -> non-sensitive).

    The shipped CSV encodes the binary outcome as +1 or -1.
    AgeInYears and the last remaining column are excluded from the design.
    A leading column of ones is prepended for the intercept.

    Parameters
    ----------
    filepath : str or Path
        Path to ``german_credit.csv``.

    Returns
    -------
    X0 : ndarray, shape (N, D)
    y0 : ndarray, shape (N,)   -- labels in {-1, +1}
    is_not_sensitive : ndarray of bool, shape (N,)
    """
    filepath = Path(filepath)
    T = pd.read_csv(filepath)

    # The first column contains the binary outcome.
    y0 = T.iloc[:, 0].values.astype(float)
    # Keep positive labels and map the other class to -1.
    y0 = np.where(y0 == 1, 1.0, -1.0)

    T = T.iloc[:, 1:]   # drop label column

    # Sensitive variable: AgeInYears > 25  -> non-sensitive
    age = T["AgeInYears"].values
    N = len(age)
    is_not_sensitive = age > 25   # bool array

    # Remove AgeInYears from predictors (sensitive attribute)
    T = T.drop(columns=["AgeInYears"])

    # Design matrix: intercept + remaining features
    X_feat = T.iloc[:, :-1].values.astype(float)  # drop last column as in original
    X0 = np.hstack([np.ones((N, 1)), X_feat])

    return X0, y0, is_not_sensitive


# ---------------------------------------------------------------------------
# Student Performance (Portuguese)
# ---------------------------------------------------------------------------

def load_student(filepath: str | Path = "student-por.csv"):
    """
    Load and preprocess the Student Performance (Portuguese) dataset.

    Sensitive variable: Sex (Male -> sensitive group, Female -> non-sensitive).
    Binary outcome: final grade G3 >= 12 -> +1 (pass), else -1 (fail).

    The feature encoding follows Code_1_Student.m exactly.

    Parameters
    ----------
    filepath : str or Path
        Path to ``student-por.csv``.

    Returns
    -------
    X0 : ndarray, shape (N, D)
    y0 : ndarray, shape (N,)  -- labels in {-1, +1}
    is_not_sensitive : ndarray of bool, shape (N,)
    """
    filepath = Path(filepath)
    T = pd.read_csv(filepath, sep=";")   # semicolon-delimited

    N = len(T)
    p = 25   # number of feature columns (without intercept)

    # ------------------------------------------------------------------
    # Map raw CSV columns to the 25-feature design used in the paper.
    # Columns are 1-indexed in the MATLAB code -> 0-indexed here.
    # ------------------------------------------------------------------
    X = np.zeros((N, p), dtype=float)

    # Numerical columns (direct copy)
    X[:, 2]  = T["age"].values
    X[:, 6]  = T["Medu"].values
    X[:, 7]  = T["Fedu"].values
    X[:, 8]  = T["traveltime"].values
    X[:, 9]  = T["studytime"].values
    X[:, 10] = T["failures"].values
    X[:, 18] = T["famrel"].values
    X[:, 19] = T["freetime"].values
    X[:, 20] = T["goout"].values
    X[:, 21] = T["Dalc"].values
    X[:, 22] = T["Walc"].values
    X[:, 23] = T["health"].values
    X[:, 24] = T["absences"].values

    # Binary-encoded categorical columns
    X[:, 0]  = (T["school"]     == "GP" ).astype(float)
    # Column index 1 (sex) is LEFT OUT as in the paper
    # (sex is the sensitive variable and is intentionally excluded)
    X[:, 3]  = (T["address"]    == "U"  ).astype(float)
    X[:, 4]  = (T["famsize"]    == "LE3").astype(float)
    X[:, 5]  = (T["Pstatus"]    == "T"  ).astype(float)
    X[:, 11] = (T["schoolsup"]  == "yes").astype(float)
    X[:, 12] = (T["famsup"]     == "yes").astype(float)
    X[:, 13] = (T["paid"]       == "yes").astype(float)
    X[:, 14] = (T["activities"] == "yes").astype(float)
    X[:, 15] = (T["nursery"]    == "yes").astype(float)
    X[:, 16] = (T["higher"]     == "yes").astype(float)
    X[:, 1]  = (T["internet"]   == "yes").astype(float)  # col 2 in MATLAB (1-indexed)
    X[:, 17] = (T["romantic"]   == "yes").astype(float)

    # Prepend intercept column
    X0 = np.hstack([np.ones((N, 1)), X])   # shape (N, 26)

    # Outcome: G3 >= 12 -> +1 (pass), else -1 (fail)
    G3 = T["G3"].values.astype(float)
    y0 = np.where(G3 >= 12, 1.0, -1.0)

    # Sensitive variable: Male -> sensitive (S=0), Female -> non-sensitive (S=1)
    is_not_sensitive = (T["sex"].values == "F")   # bool array

    return X0, y0, is_not_sensitive


def load_acs_income(
    state: str = "NV",
    year:  str = "2018",
    data_dir: str = "data",
    income_threshold: float = 50_000.0,
):
    try:
        from folktables import ACSDataSource, ACSIncome
    except ImportError:
        raise ImportError("Install folktables:  pip install folktables")

    data_source = ACSDataSource(
        survey_year = year,
        horizon     = "1-Year",
        survey      = "person",
        root_dir    = str(data_dir),
    )
    acs_data = data_source.get_data(states=[state], download=True)

    X_raw, y_raw, _ = ACSIncome.df_to_numpy(acs_data)

    y0 = np.where(y_raw, 1.0, -1.0)

    sex_col          = list(ACSIncome.features).index("SEX")
    sex              = X_raw[:, sex_col]
    is_not_sensitive = (sex == 2)

    X_feat = np.nan_to_num(X_raw, nan=-1.0).astype(float)
    N      = X_feat.shape[0]
    X0     = np.hstack([np.ones((N, 1)), X_feat])

    print(f"ACSIncome {state} {year}: N={N}  D={X0.shape[1]}  "
          f"pos_rate={y0.mean()*100:.1f}%  "
          f"sensitive(Male)={(~is_not_sensitive).mean()*100:.1f}%")

    return X0, y0, is_not_sensitive
