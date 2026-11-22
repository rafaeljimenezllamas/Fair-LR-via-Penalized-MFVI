"""
main.py
-------
Main entry point for the Fair MFVI logistic regression pipeline.

Mirrors ``Code_5.m``.

Workflow
--------
1. Load and preprocess a dataset (German Credit or Student Performance).
2. Shuffle the data with a fixed random seed for reproducibility.
3. For each of the K = 5 cross-validation folds:
   a. Create train / test split and compute the fairness vectors z, z_T.
   b. Search hyperparameters (a0, b0) for each lambda (``optimizer.run_fold``).
   c. Save per-fold results to runs/Results_<DATASET>_CV_<n>.npz
4. Print total runtime.

Usage
-----
    python main.py --dataset german
    python main.py --dataset student
    python main.py --dataset german --num_folds 5 --enforce_on test
"""

import argparse
import math
import time
from pathlib import Path

import numpy as np

from data_loader       import load_german, load_student
from cross_validation  import train_test_split_cv
from optimizer         import run_fold


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fair MFVI logistic regression – main pipeline"
    )
    parser.add_argument(
        "--dataset",
        choices=["german", "student"],
        default="german",
        help="Dataset to use (default: german).",
    )
    parser.add_argument(
        "--num_folds",
        type=int,
        default=5,
        help="Number of cross-validation folds (default: 5).",
    )
    parser.add_argument(
        "--enforce_on",
        choices=["train", "test"],
        default="test",
        help=(
            "Whether to enforce the fairness constraint using the train or "
            "test set to compute z (default: test)."
        ),
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "data"),
        help="Directory containing the input CSV files (default: repository data/).",
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "runs"),
        help="Directory for new per-fold results (default: repository runs/).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for data shuffling (default: 42).",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Load dataset
    # ------------------------------------------------------------------ #
    data_dir = Path(args.data_dir)
    if args.dataset == "german":
        X0, y0, is_not_sensitive = load_german(data_dir / "german_credit.csv")
        tag = "GERMAN"
    else:
        X0, y0, is_not_sensitive = load_student(data_dir / "student-por.csv")
        tag = "STUDENT"

    N = len(y0)
    print(f"Dataset: {tag}  |  N={N}  |  D={X0.shape[1]}")

    # Check that the design matrix is full rank
    if np.linalg.matrix_rank(X0) < X0.shape[1]:
        raise ValueError("Design matrix X0 does not have full column rank.")

    # ------------------------------------------------------------------ #
    # 2. Shuffle data (fixed seed for reproducibility)
    # ------------------------------------------------------------------ #
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(N)
    X0              = X0[perm]
    y0              = y0[perm]
    is_not_sensitive = is_not_sensitive[perm]

    # ------------------------------------------------------------------ #
    # 3. Define lambda grid (60 log-spaced values from 1e-3 to 1e8)
    # ------------------------------------------------------------------ #
    vec_lambda = np.logspace(-3, 8, 60)

    # Approximate test-set size per fold
    M = math.floor(N / args.num_folds)

    # ------------------------------------------------------------------ #
    # 4. Cross-validation loop
    # ------------------------------------------------------------------ #
    t0_total = time.time()

    for n in range(1, args.num_folds + 1):
        print(f"\n{'='*60}")
        print(f"  Fold {n} / {args.num_folds}")
        print(f"{'='*60}")

        # Create train / test split and fairness vectors
        (X, y, z, XT, yT, zT,
         ind_train, ind_test, *_) = train_test_split_cv(
            X0, y0, is_not_sensitive, M, n
        )

        # Determine fairness direction (train z or test zT)
        enforce_on = args.enforce_on   # "train" or "test"

        # Output file for this fold
        save_path = results_dir / f"Results_{tag}_CV_{n}.npz"

        t0_fold = time.time()
        run_fold(
            X=X, y=y, z=z,
            XT=XT, yT=yT, zT=zT,
            enforce_fairness_on=enforce_on,
            vec_lambda=vec_lambda,
            save_path=save_path,
        )
        elapsed_fold = time.time() - t0_fold
        print(f"  Fold {n} done in {elapsed_fold/60:.1f} min")

    elapsed_total = time.time() - t0_total
    print(f"\nAll {args.num_folds} folds completed in {elapsed_total/60:.1f} min")
    print(f"Results saved in: {results_dir.resolve()}")


if __name__ == "__main__":
    main()
