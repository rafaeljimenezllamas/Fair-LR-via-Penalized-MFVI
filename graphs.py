"""Plot saved fair MFVI results without fitting any models."""

import argparse
import math
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from cross_validation import train_test_split_cv
from data_loader import load_german
from h_beta_hist import h_beta_histograms
from load_acs_pvar import (
    build_design_matrix,
    load_acs_income_full,
    total_dummy_columns,
)
from utils import logdet


ROOT = Path(__file__).resolve().parent
SAVED = ROOT / "results"
ACS_TAG = "ACS_NV_P727_0"
QUICK_FIGURES = ("1", "2", "3", "4", "5-right", "6", "7", "S1", "S2", "S4")
CURVE_FIELDS = (
    "sr_train", "sr_test", "p2_5_train", "p97_train",
    "p2_5_test", "p97_test", "dif_test",
)


def save(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(directory / f"{name}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


def fold_paths(directory, tag):
    return [directory / f"Results_{tag}_CV_{fold}.npz" for fold in range(1, 6)]


def load_curves(directory, tag):
    values = {key: [] for key in CURVE_FIELDS}
    grid = None
    for path in fold_paths(directory, tag):
        with np.load(path, allow_pickle=False) as fold:
            lam = fold["vec_lambda"]
            if grid is None:
                grid = lam
            if not np.array_equal(lam, grid):
                raise ValueError(f"Penalty grid differs between folds in {path}")
            for key in values:
                value = fold[key]
                if value.shape != lam.shape:
                    raise ValueError(f"{key} has the wrong shape in {path}")
                values[key].append(value)
    return grid, {key: np.column_stack(val).mean(axis=1) for key, val in values.items()}


def success_rate(ax, lam, curves, intervals=True):
    for split, label, color in (("train", "Train", "r"), ("test", "Test", "b")):
        ax.semilogx(lam, curves[f"sr_{split}"], color=color, label=f"{label} Set")
        if intervals:
            ax.semilogx(lam, curves[f"p2_5_{split}"], "--", color=color,
                        label=f"95% PI ({label})")
            ax.semilogx(lam, curves[f"p97_{split}"], "--", color=color)
    ax.set(xlabel=r"$\lambda$", ylabel="Success Rate")
    ax.legend()


def test_unfairness(ax, lam, curves):
    ax.loglog(lam, curves["dif_test"], "--b", label="Test Set")
    ax.set(xlabel=r"$\lambda$", ylabel=r"$\mathcal{F}$")
    ax.legend()


def plot_tradeoff(directory, tag, out, supplement=False):
    lam, curves = load_curves(directory, tag)
    panels = [(success_rate, "left"), (test_unfairness, "right")]
    if supplement:
        panels = [
            (success_rate, "top-left"),
            (lambda ax, x, y: success_rate(ax, x, y, intervals=False), "top-right"),
            (test_unfairness, "bottom"),
        ]
    number = {"GERMAN": "03", ACS_TAG: "04", "STUDENT": "S1"}[tag]
    for draw, position in panels:
        fig, ax = plt.subplots(figsize=(8, 5))
        draw(ax, lam, curves)
        ax.grid(True, which="both", ls=":")
        save(fig, out, f"figure-{number}-{position}")


def fold_objective(path, design):
    """Re-evaluate the notebook ELBO using the first saved Gamma parameters."""
    X, y, z, _, y_test, z_test, *_ = design
    # The original posterior histories store NumPy object arrays.
    with np.load(path, allow_pickle=True) as fold:
        if not (np.array_equal(fold["y"], y) and np.array_equal(fold["yT"], y_test)
                and np.allclose(fold["z"], z) and np.allclose(fold["zT"], z_test)):
            raise ValueError(f"Raw data and saved split differ in {path}")
        betas, sigmas = fold["best_beta"], fold["best_Sigma"]
        va, vb = fold["best_vec_a"], fold["best_vec_b"]

    a0, b0 = float(va[0][0]), float(vb[0][0])
    ak, bk = float(va[0][-1]), float(vb[0][-1])
    objective = np.zeros(len(betas))
    for i, entry in enumerate(betas):
        if entry is None:
            continue
        beta = entry.astype(float)
        sigma = sigmas[i].astype(float)
        inv_sigma = np.linalg.inv(sigma)
        xsx = (X @ sigma * X).sum(axis=1)
        xi = np.sqrt(np.maximum(xsx + (X @ beta) ** 2, 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            omega = np.tanh(xi / 2.0) / (4.0 * xi)
        omega = np.where(np.isnan(omega), 1.0 / 8.0, omega)
        lchi = -np.log1p(np.exp(-xi)) - xi / 2.0 + xi**2 * omega
        objective[i] = (
            0.5 * (beta @ inv_sigma @ beta + logdet(sigma))
            + a0 * np.log(b0) - ak * np.log(bk)
            + math.lgamma(ak) - math.lgamma(a0)
            + ak - ak * b0 / bk + lchi.sum()
        )
    return objective


def mean_objective(directory, tag, data):
    X, y, group = data
    order = np.random.default_rng(42).permutation(len(y))
    X, y, group = X[order], y[order], group[order]
    paths = fold_paths(directory, tag)
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    return np.column_stack([
        fold_objective(path, train_test_split_cv(X, y, group, len(y) // 5, fold))
        for fold, path in enumerate(paths, 1)
    ]).mean(axis=1)


def plot_pareto(out, number, curves, objective=None):
    if objective is None:
        fig, ax = plt.subplots(figsize=(6, 5))
        panels = [(ax, curves["sr_test"], "Estimated Success Rate")]
    else:
        fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5))
        panels = [
            (left, objective, r"$\mathcal{L}^*(q^*,\zeta^*)$"),
            (right, curves["sr_test"], "Estimated Success Rate"),
        ]
    fairness = curves["dif_test"]
    valid = fairness > 0
    if objective is not None:
        valid &= np.isfinite(objective)
    order = np.argsort(fairness[valid])
    for ax, values, label in panels:
        ax.semilogx(fairness[valid][order], values[valid][order], "k")
        ax.set(xlabel=r"$\mathcal{F}$", ylabel=label)
        ax.grid(True, which="both", ls=":")
    name = f"figure-{number}" if objective is not None else f"figure-{number}-right"
    save(fig, out, name)


def plot_surface(out, saved):
    with np.load(saved / "GERMAN_hyperparameter_surface.npz", allow_pickle=False) as data:
        x, y = np.meshgrid(np.log10(data["a0_grid"]), np.log10(data["ratio_grid"]),
                           indexing="ij")
        sr = data["SR"]
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": "3d"})
    surface = ax.plot_surface(x, y, sr, cmap="viridis")
    ticks = [-5, 0, 5, 10, 15]
    ax.set(xlabel=r"$a_0$", ylabel=r"$a_0/b_0$", title="Success Rate (Train)",
           xticks=ticks, yticks=ticks)
    ax.view_init(elev=26, azim=45)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda value, _: rf"$10^{{{value:g}}}$"))
    fig.colorbar(surface, ax=ax, shrink=0.65, pad=0.10)
    save(fig, out, "figure-01")


def plot_runtime(out, saved):
    rows = pd.read_csv(saved / "ACS_NV_P727_0_runtime_vs_p_lambda.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    for lam, points in rows.groupby("lambda", sort=True):
        points = points.sort_values("p")
        ax.loglog(points["p"], points["time_median"], "o-", label=rf"$\lambda={lam:.1e}$")
    ax.set(xlabel="Design columns (including intercept)", ylabel="Median CAVI runtime (s)")
    ax.legend()
    ax.grid(True, which="both", ls=":")
    save(fig, out, "figure-06")


def plot_comparison(out, saved, dataset, number):
    with np.load(saved / "comparisons" / f"{dataset}.npz", allow_pickle=False) as data:
        lam = data["vec_lambda"]
        series = {
            key: data[key].mean(axis=0) for key in (
                "bayesian_accuracy", "bayesian_unfairness",
                "frequentist_accuracy", "frequentist_unfairness",
            )
        }
    if len(lam) != 61 or lam[0] != 0 or np.any(lam[1:] <= 0):
        raise ValueError(f"Expected one zero penalty and 60 positive penalties for {dataset}")

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.5))
    palette = (("frequentist", "Frequentist", "#009E73"),
               ("bayesian", "Bayesian", "#E69F00")) if dataset == "german" else (
                   ("frequentist", "Frequentist", "#2166AC"),
                   ("bayesian", "Bayesian", "#B2182B"),
               )
    for mode, label, color in palette:
        accuracy = series[f"{mode}_accuracy"]
        unfairness = series[f"{mode}_unfairness"]
        left.axhline(accuracy[0], ls="--", color=color, label=f"Standard {label} LR")
        left.semilogx(lam[1:], accuracy[1:], color=color, label=f"Fair {label} LR")
        right.axhline(unfairness[0], ls="--", color=color, label=f"Standard {label} LR")
        right.loglog(lam[1:], unfairness[1:], color=color, label=f"Fair {label} LR")
    left.set(ylabel="Test accuracy")
    right.set(ylabel="Unfairness")
    for ax in (left, right):
        ax.set(xlabel=r"$\lambda$")
        ax.legend(fontsize="small")
        ax.grid(True, which="both", ls=":")
    save(fig, out, f"figure-{number}")


def _load_fold(results_dir: Path, tag: str, n: int) -> dict:
    """Load a single fold .npz file and return its contents as a dict."""
    path = results_dir / f"Results_{tag}_CV_{n}.npz"
    with np.load(path, allow_pickle=True) as data:
        return dict(data)


def plot_fold_diagnostics(args):
    """Plot the original train/test and posterior diagnostics for saved folds."""
    tag         = args.dataset.upper()
    results_dir = Path(args.results_dir)
    out_dir     = Path(args.out_dir) / tag.lower()
    n_folds     = args.num_folds

    # ------------------------------------------------------------------ #
    # Aggregate results across folds
    # ------------------------------------------------------------------ #
    fold0 = _load_fold(results_dir, tag, 1)
    vec_lambda = fold0["vec_lambda"]
    n_lam      = len(vec_lambda)

    sr_train_cv   = np.zeros((n_lam, n_folds))
    sr_test_cv    = np.zeros((n_lam, n_folds))
    p97_train_cv  = np.zeros((n_lam, n_folds))
    p2_5_train_cv = np.zeros((n_lam, n_folds))
    p97_test_cv   = np.zeros((n_lam, n_folds))
    p2_5_test_cv  = np.zeros((n_lam, n_folds))
    dif_train_cv  = np.zeros((n_lam, n_folds))
    dif_test_cv   = np.zeros((n_lam, n_folds))

    # Per-fold storage for H(beta) histograms
    H_beta_train_folds = []
    H_beta_test_folds  = []

    for n in range(1, n_folds + 1):
        fold = _load_fold(results_dir, tag, n)

        sr_train_cv[:,   n-1] = fold["sr_train"]
        sr_test_cv[:,    n-1] = fold["sr_test"]
        p97_train_cv[:,  n-1] = fold["p97_train"]
        p2_5_train_cv[:, n-1] = fold["p2_5_train"]
        p97_test_cv[:,   n-1] = fold["p97_test"]
        p2_5_test_cv[:,  n-1] = fold["p2_5_test"]
        dif_train_cv[:,  n-1] = fold["dif_train"]
        dif_test_cv[:,   n-1] = fold["dif_test"]

        # -- Compute H(beta) distribution via Monte Carlo (Eq. 20) --
        #    H(beta) = (beta^T z)^2 = beta^T (z z^T) beta
        #    We sample beta ~ N(beta*, Sigma*) and evaluate H for each lambda.
        n_sim  = 10_000
        betas  = fold["best_beta"]   # array of (D,) arrays, one per lambda
        sigmas = fold["best_Sigma"]  # array of (D,D) arrays
        z      = fold["z"]
        zT     = fold["zT"]

        H_train = np.zeros((n_sim, n_lam))
        H_test  = np.zeros((n_sim, n_lam))

        for i in range(n_lam):
            if betas[i] is None:
                continue
            b  = np.asarray(betas[i], dtype=float)
            S  = np.asarray(sigmas[i], dtype=float)
            # Sample from posterior: shape (n_sim, D)
            samples = np.random.multivariate_normal(b, S, size=n_sim)
            # H(beta) = (beta^T z)^2   (scalar per sample)
            H_train[:, i] = (samples @ z)  ** 2
            H_test[:,  i] = (samples @ zT) ** 2

        H_beta_train_folds.append(H_train)
        H_beta_test_folds.append(H_test)

    # Cross-validated means
    sr_train_mean   = sr_train_cv.mean(axis=1)
    sr_test_mean    = sr_test_cv.mean(axis=1)
    p97_train_mean  = p97_train_cv.mean(axis=1)
    p2_5_train_mean = p2_5_train_cv.mean(axis=1)
    p97_test_mean   = p97_test_cv.mean(axis=1)
    p2_5_test_mean  = p2_5_test_cv.mean(axis=1)
    dif_train_mean  = dif_train_cv.mean(axis=1)
    dif_test_mean   = dif_test_cv.mean(axis=1)

    fs   = 14    # base font size
    lw   = 2     # line width

    # ------------------------------------------------------------------ #
    # Figure 1: Success Rate vs lambda (with 95% PI)
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogx(vec_lambda, sr_train_mean,   "r",   lw=lw,  label="Train Set")
    ax.semilogx(vec_lambda, sr_test_mean,    "b",   lw=lw,  label="Test Set")
    ax.semilogx(vec_lambda, p97_train_mean,  "--r", lw=lw,  label="95% PI (Train)")
    ax.semilogx(vec_lambda, p2_5_train_mean, "--r", lw=lw)
    ax.semilogx(vec_lambda, p97_test_mean,   "--b", lw=lw,  label="95% PI (Test)")
    ax.semilogx(vec_lambda, p2_5_test_mean,  "--b", lw=lw)
    ax.set_xlabel(r"$\lambda$", fontsize=24)
    ax.set_ylabel("Success Rate", fontsize=24)
    ax.tick_params(labelsize=fs)
    ax.legend(fontsize=16)
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    save(fig, out_dir, f"{tag}_fig1_success_rate_with_PI")

    # ------------------------------------------------------------------ #
    # Figure 2: Success Rate vs lambda (without PI)
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogx(vec_lambda, sr_train_mean, "r", lw=lw, label="Train Set")
    ax.semilogx(vec_lambda, sr_test_mean,  "b", lw=lw, label="Test Set")
    ax.set_xlabel(r"$\lambda$", fontsize=24)
    ax.set_ylabel("Success Rate", fontsize=24)
    ax.tick_params(labelsize=fs)
    ax.legend(fontsize=16)
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    save(fig, out_dir, f"{tag}_fig2_success_rate")

    # ------------------------------------------------------------------ #
    # Figure 3: Unfairness metric F vs lambda (log-log)
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(vec_lambda, dif_train_mean, "--r", lw=lw, label="Train Set")
    ax.loglog(vec_lambda, dif_test_mean,  "--b", lw=lw, label="Test Set")
    ax.set_xlabel(r"$\lambda$", fontsize=24)
    ax.set_ylabel(r"$\mathcal{F}_0$", fontsize=24)
    ax.tick_params(labelsize=fs)
    ax.legend(fontsize=16)
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    save(fig, out_dir, f"{tag}_fig3_unfairness_metric")

    # ------------------------------------------------------------------ #
    # Figures 4-6: H(beta) histograms for lambda indices 0, 32, 59
    # (approximately 1e-3, middle, 1e8)
    # Use the last fold (fold 5 in MATLAB, index n_folds-1 here)
    # ------------------------------------------------------------------ #
    H_train_last = H_beta_train_folds[-1]   # (n_sim, n_lam)
    H_test_last  = H_beta_test_folds[-1]

    for fig_num, lam_idx in enumerate([0, 32, 59], start=4):
        lam_val = vec_lambda[lam_idx]
        h_tr    = H_train_last[:, lam_idx]
        h_te    = H_test_last[:,  lam_idx]

        fig, ax = plt.subplots(figsize=(7, 4))
        # Use log-spaced bins for the x-axis as in the paper
        all_vals = np.concatenate([h_tr, h_te])
        all_vals = all_vals[all_vals > 0]
        bins = np.logspace(np.log10(all_vals.min()), np.log10(all_vals.max()), 40)

        ax.hist(h_tr, bins=bins, density=True, alpha=0.6, color="r", label="Train Set")
        ax.hist(h_te, bins=bins, density=True, alpha=0.6, color="b", label="Test Set")
        ax.set_xscale("log")
        ax.set_xlabel(r"$H(\beta)$", fontsize=24)
        ax.set_ylabel(r"$f(\cdot)$", fontsize=24)
        ax.set_title(rf"$\lambda = {lam_val:.1e}$", fontsize=fs)
        ax.tick_params(labelsize=fs)
        ax.legend(fontsize=16)
        ax.grid(True, which="both", ls=":")
        fig.tight_layout()
        save(fig, out_dir, f"{tag}_fig{fig_num}_H_beta_lam{lam_idx}")

    # ------------------------------------------------------------------ #
    # Figure 7: Posterior mean coefficients vs lambda
    # (first fold, first 6 non-intercept coefficients)
    # ------------------------------------------------------------------ #
    fold1 = _load_fold(results_dir, tag, 1)
    betas = fold1["best_beta"]   # (n_lam,) array of (D,) arrays

    n_coefs = min(6, betas[0].shape[0] - 1)   # exclude intercept
    beta_traj = np.zeros((n_lam, n_coefs))
    for i in range(n_lam):
        if betas[i] is not None:
            beta_traj[i] = betas[i][1 : 1 + n_coefs]

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, n_coefs))
    for j in range(n_coefs):
        ax.semilogx(vec_lambda, beta_traj[:, j], lw=lw, color=colors[j],
                    label=rf"$\beta_{{{j+1}}}$")
    ax.set_xlabel(r"$\lambda$", fontsize=24)
    ax.set_ylabel("Values", fontsize=24)
    ax.tick_params(labelsize=fs)
    ax.legend(fontsize=14, ncol=2)
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    save(fig, out_dir, f"{tag}_fig7_beta_evolution")

    print(f"\nAll figures saved to: {out_dir.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--figure", choices=("all", *QUICK_FIGURES, "5"), default="all",
                      help="Use 5 for full ACS ELBO replay with the complete fold archives.")
    mode.add_argument("--dataset", choices=("german", "student"),
                      help="Draw diagnostics from a directory of saved experiment folds.")
    parser.add_argument("--results-dir", "--results_dir", type=Path, default=ROOT / "runs")
    parser.add_argument("--num-folds", "--num_folds", type=int, default=5)
    parser.add_argument("--saved-dir", type=Path, default=SAVED)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--out-dir", "--out_dir", type=Path, default=ROOT / "runs" / "figures")
    parser.add_argument("--acs-fold-dir", type=Path)
    parser.add_argument("--acs-data-dir", type=Path, default=ROOT / "data",
                        help="Folktables root that already contains the verified Nevada CSV.")
    args = parser.parse_args()
    if args.out_dir.resolve().is_relative_to((ROOT / "figures").resolve()):
        parser.error("Choose an output directory outside the published figures/ directory")
    if args.dataset:
        plot_fold_diagnostics(args)
        return
    if args.figure == "5":
        if args.acs_fold_dir is None:
            parser.error("Figure 5 requires --acs-fold-dir")
        if not (args.acs_data_dir / "2018" / "1-Year" / "psam_p32.csv").is_file():
            parser.error("--acs-data-dir must contain 2018/1-Year/psam_p32.csv")
        missing = [p for p in fold_paths(args.acs_fold_dir, ACS_TAG) if not p.is_file()]
        if missing:
            parser.error(f"Missing ACS fold archive: {missing[0]}")
    figures = QUICK_FIGURES if args.figure == "all" else (args.figure,)
    saved, out = args.saved_dir, args.out_dir
    for number in figures:
        print(f"Redrawing Figure {number}", flush=True)
        if number == "1":
            plot_surface(out, saved)
        elif number == "2":
            directory = saved / "folds"
            _, curves = load_curves(directory, "GERMAN")
            data = load_german(args.data_dir / "german_credit.csv")
            plot_pareto(out, "02", curves, mean_objective(directory, "GERMAN", data))
        elif number == "3":
            plot_tradeoff(saved / "folds", "GERMAN", out)
        elif number == "4":
            plot_tradeoff(saved / "acs_curves", ACS_TAG, out)
        elif number in ("5", "5-right"):
            directory = args.acs_fold_dir if number == "5" else saved / "acs_curves"
            _, curves = load_curves(directory, ACS_TAG)
            objective = None
            if number == "5":
                numeric, dummies, y, group = load_acs_income_full(data_dir=args.acs_data_dir)
                design = build_design_matrix(numeric, dummies, total_dummy_columns(dummies))
                objective = mean_objective(args.acs_fold_dir, ACS_TAG, (design, y, group))
            plot_pareto(out, "05", curves, objective)
        elif number == "6":
            plot_runtime(out, saved)
        elif number == "7":
            plot_comparison(out, saved, "german", "07")
        elif number == "S1":
            plot_tradeoff(saved / "folds", "STUDENT", out, supplement=True)
        elif number == "S2":
            plot_comparison(out, saved, "student", "S2")
        elif number == "S4":
            for tag in ("GERMAN", "STUDENT"):
                fig = h_beta_histograms(
                    tag, lam_indices=(0, 32, 59), fold=1, which="test",
                    n_sim=20_000, seed=0, results_dir=saved / "s4",
                    figures_dir=out / "figure-S4", show=False,
                )
                plt.close(fig)


if __name__ == "__main__":
    main()
