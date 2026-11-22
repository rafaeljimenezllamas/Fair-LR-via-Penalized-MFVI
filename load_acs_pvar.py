"""
load_acs_pvar.py
-----------------
Load ACSIncome with full one-hot encoding and build design matrices
from prefixes of the ordered dummy columns.

SEX defines the sensitive group and is excluded from the design.
"""

import numpy as np
import pandas as pd


NUMERIC_COLS = ["AGEP", "WKHP"]
# Dummy blocks are ordered by increasing cardinality.
CATEGORICAL_COLS_ORDER = ["COW", "MAR", "RELP", "RAC1P", "SCHL", "POBP", "OCCP"]


def load_acs_income_full(state: str = "NV", year: str = "2018", data_dir: str = "data"):
    """
    Load ACSIncome and return the numeric columns (AGEP, WKHP), dummy
    blocks in CATEGORICAL_COLS_ORDER, labels in {+1, -1}, and the
    non-sensitive-group mask (SEX == 2).
    """
    from folktables import ACSDataSource, ACSIncome

    data_source = ACSDataSource(
        survey_year=year, horizon="1-Year", survey="person", root_dir=str(data_dir)
    )
    acs_data = data_source.get_data(states=[state], download=True)

    df = ACSIncome._preprocess(acs_data)
    y_raw = ACSIncome.target_transform(df[ACSIncome.target]).to_numpy()
    y0 = np.where(y_raw, 1.0, -1.0)

    sex = df["SEX"].to_numpy()
    is_not_sensitive = (sex == 2)

    df_numeric = df[NUMERIC_COLS].astype(float).reset_index(drop=True)

    dummy_blocks = {}
    for col in CATEGORICAL_COLS_ORDER:
        d = pd.get_dummies(
            df[col].astype("category"), prefix=col, drop_first=False
        ).astype(float).reset_index(drop=True)
        dummy_blocks[col] = d

    return df_numeric, dummy_blocks, y0, is_not_sensitive


def build_design_matrix(df_numeric, dummy_blocks, n_dummy_cols: int):
    """
    Build X0 from an intercept, the numeric columns, and the first
    n_dummy_cols columns in the concatenated, ordered dummy blocks.
    """
    all_dummies = pd.concat([dummy_blocks[c] for c in CATEGORICAL_COLS_ORDER], axis=1)
    n_dummy_cols = min(n_dummy_cols, all_dummies.shape[1])
    dummies_subset = all_dummies.iloc[:, :n_dummy_cols]

    feat = pd.concat([df_numeric, dummies_subset], axis=1)
    X_feat = feat.to_numpy(dtype=float)
    N = X_feat.shape[0]
    X0 = np.hstack([np.ones((N, 1)), X_feat])
    return X0


def total_dummy_columns(dummy_blocks):
    return sum(d.shape[1] for d in dummy_blocks.values())
