"""
comparison.py
-------------
Standalone evaluation module for the Fair MFVI paper revision.

Implements the original three-method evaluation and supporting metrics.
The four-method paper figures use the archived comparison arrays.

Usage from notebook:
    import comparison, importlib; importlib.reload(comparison)
    results = comparison.run_comparison(dataset='german', ...)
"""

import math
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, f1_score,
    balanced_accuracy_score, brier_score_loss,
)
from sklearn.preprocessing import StandardScaler

from data_loader import load_german, load_student
from cross_validation import train_test_split_cv


# ---------------------------------------------------------------------------
# Fairness metrics (standard, for reporting to reviewer)
# ---------------------------------------------------------------------------

def demographic_parity(y_pred_binary, sensitive):
    p_s  = y_pred_binary[sensitive].mean()
    p_ns = y_pred_binary[~sensitive].mean()
    return float(abs(p_ns - p_s))


def disparate_impact(y_pred_binary, sensitive):
    p_s  = y_pred_binary[sensitive].mean()
    p_ns = y_pred_binary[~sensitive].mean()
    if p_ns < 1e-9:
        return np.nan
    return float(p_s / p_ns)


def equal_opportunity(y_pred_binary, y_true, sensitive):
    def tpr(mask):
        pos = y_true[mask] == 1
        if pos.sum() == 0:
            return np.nan
        return float(y_pred_binary[mask][pos].mean())
    t_s, t_ns = tpr(sensitive), tpr(~sensitive)
    if np.isnan(t_s) or np.isnan(t_ns):
        return np.nan
    return float(abs(t_ns - t_s))


def compute_all_metrics(y_true, y_pred_binary, y_pred_proba, sensitive):
    yt01 = (y_true        == 1).astype(int)
    yp01 = (y_pred_binary == 1).astype(int)
    return dict(
        accuracy     = float((y_pred_binary == y_true).mean()),
        balanced_acc = balanced_accuracy_score(yt01, yp01),
        auc          = roc_auc_score(yt01, y_pred_proba),
        f1           = f1_score(yt01, yp01, zero_division=0),
        brier        = brier_score_loss(yt01, y_pred_proba),
        dem_parity   = demographic_parity(yp01, sensitive),
        disp_impact  = disparate_impact(yp01, sensitive),
        eq_opp       = equal_opportunity(yp01, yt01, sensitive),
    )


# ---------------------------------------------------------------------------
# Fairness proxy for each method
# ---------------------------------------------------------------------------

def fairness_proxy_freq(beta, z):
    """
    Frequentist fairness proxy: (beta^T z)^2.
    Point estimate — no uncertainty over beta.
    """
    return float((beta @ z) ** 2)


def fairness_proxy_mfvi(beta, Sigma, z):
    """
    MFVI fairness metric (Eq. 29): (beta*^T z)^2 + Tr(z^T Sigma* z).
    Includes the uncertainty contribution from the posterior covariance.
    """
    return float((beta @ z) ** 2 + float(z @ Sigma @ z))


# ---------------------------------------------------------------------------
# Safe cast of object arrays loaded from .npz
# ---------------------------------------------------------------------------

def _to_float(arr):
    if arr is None:
        return None
    try:
        return np.array(arr, dtype=np.float64)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Standard LR (no fairness)
# ---------------------------------------------------------------------------

def _fit_standard_lr(X_tr, y_tr, X_te):
    scaler = StandardScaler()
    Xtr_s  = scaler.fit_transform(X_tr[:, 1:])
    Xte_s  = scaler.transform(X_te[:, 1:])
    clf = LogisticRegression(max_iter=1000, solver="lbfgs", C=1.0)
    clf.fit(Xtr_s, y_tr)
    proba = clf.predict_proba(Xte_s)[:, 1]
    preds = np.where(proba >= 0.5, 1.0, -1.0)
    return preds, proba


# ---------------------------------------------------------------------------
# Core CV loop
# ---------------------------------------------------------------------------

def run_comparison(
    dataset      = "german",
    data_dir     = ".",
    results_dir  = "results",
    figures_dir  = "figures",
    num_folds    = 5,
    seed         = 42,
):
    """
    Run the full comparison experiment with lambda sweep for all methods.

    For each fold and each lambda:
      - Standard LR: single fit (no lambda dependence), shown as horizontal line
      - Fair LR (freq): penalized with (lam/2)*(beta^T z)^2, lambda sweep
      - MFVI (ours): loaded from saved .npz, lambda sweep

    Fairness proxies reported:
      - Fair LR (freq): (beta^T z)^2            [no posterior uncertainty]
      - MFVI (ours):    (beta*^T z)^2 + Tr(z^T Sigma* z)  [Eq. 29]

    These are NOT the same quantity — clearly noted in plots and tables.
    """
    from vb_logit_pred import vb_logit_pred

    data_dir    = Path(data_dir)
    results_dir = Path(results_dir)
    tag         = dataset.upper()
    out_dir     = Path(figures_dir) / dataset.lower()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  COMPARISON: {tag}")
    print(f"{'='*60}")

    # ------------------------------------------------------------------ #
    # Load data — same shuffle as original pipeline
    # ------------------------------------------------------------------ #

    if dataset == "german":
        X0, y0, is_not_sensitive = load_german(data_dir / "german_credit.csv")
    elif dataset == "student":
        X0, y0, is_not_sensitive = load_student(data_dir / "student-por.csv")
    else:
        # ACS dataset: extract state from tag e.g. 'acs_NV' -> 'NV'
        state = dataset.split("_")[1].upper()
        from data_loader import load_acs_income
        X0, y0, is_not_sensitive = load_acs_income(state=state, data_dir=str(data_dir))

    N = len(y0)
    print(f"  N={N}  D={X0.shape[1]}")

    rng  = np.random.default_rng(seed)
    perm = rng.permutation(N)
    X0               = X0[perm]
    y0               = y0[perm]
    is_not_sensitive = is_not_sensitive[perm]
    M = math.floor(N / num_folds)

    # ------------------------------------------------------------------ #
    # Lambda grid (same as original pipeline)
    # ------------------------------------------------------------------ #
    fold0_data = np.load(results_dir / f"Results_{tag}_CV_1.npz",
                         allow_pickle=True)
    vec_lambda = fold0_data["vec_lambda"]
    n_lam      = len(vec_lambda)
    print(f"  Lambda grid: {n_lam} values  "
          f"[{vec_lambda[0]:.1e} .. {vec_lambda[-1]:.1e}]")

    metric_keys = ["accuracy", "balanced_acc", "auc", "f1", "brier",
                   "dem_parity", "disp_impact", "eq_opp"]

    # Storage — shape (n_folds, n_lam)
    acc_lr_folds    = np.full((num_folds, n_lam), np.nan)  # constant per fold
    acc_freq_folds  = np.full((num_folds, n_lam), np.nan)
    acc_mfvi_folds  = np.full((num_folds, n_lam), np.nan)

    proxy_freq_folds = np.full((num_folds, n_lam), np.nan)
    proxy_mfvi_folds = np.full((num_folds, n_lam), np.nan)

    # Extended metrics at each lambda for MFVI — (n_folds, n_lam, n_metrics)
    mfvi_ext_folds     = np.full((num_folds, n_lam, len(metric_keys)), np.nan)
    mfvi_ext_tr_folds  = np.full((num_folds, n_lam, len(metric_keys)), np.nan)

    # Runtime (seconds per fold)
    rt_lr_folds   = np.zeros(num_folds)
    rt_freq_folds = np.zeros(num_folds)
    rt_mfvi_folds = np.zeros(num_folds)

    for fold_n in range(1, num_folds + 1):
        print(f"\n  -- Fold {fold_n}/{num_folds} --")

        (X_tr, y_tr, z_tr,
         X_te, y_te, z_te,
         idx_tr, idx_te, *_) = train_test_split_cv(
            X0, y0, is_not_sensitive, M, fold_n
        )
        sens_te = ~is_not_sensitive[idx_te]
        sens_tr = ~is_not_sensitive[idx_tr]

        fold_data = np.load(
            results_dir / f"Results_{tag}_CV_{fold_n}.npz",
            allow_pickle=True
        )

        # ---- Standard LR (no fairness, single fit) ----
        t0 = time.time()
        preds_lr, proba_lr = _fit_standard_lr(X_tr, y_tr, X_te)
        rt_lr_folds[fold_n-1] = time.time() - t0
        acc_lr = float((preds_lr == y_te).mean())
        acc_lr_folds[fold_n-1, :] = acc_lr   # same value for all lambdas
        print(f"    Standard LR:  acc={acc_lr:.3f}  t={rt_lr_folds[fold_n-1]:.2f}s")

        # ---- Fair LR (freq) — lambda sweep ----
        print(f"    Fair LR (freq): sweeping {n_lam} lambdas...", end=" ", flush=True)
        scaler = StandardScaler()
        Xtr_s  = scaler.fit_transform(X_tr[:, 1:])
        Xte_s  = scaler.transform(X_te[:, 1:])
        N_tr, N_te = Xtr_s.shape[0], Xte_s.shape[0]
        Xtr_s = np.hstack([np.ones((N_tr, 1)), Xtr_s])
        Xte_s = np.hstack([np.ones((N_te, 1)), Xte_s])
        D     = Xtr_s.shape[1]
        z_dir = z_te.ravel()   # use test z (same convention as MFVI)

        t0 = time.time()
        for i, lam in enumerate(vec_lambda):
            beta = np.zeros(D)
            lr   = min(0.01, 1.0 / (1.0 + lam))   # smaller lr for large lam
            for _ in range(500):
                margins = y_tr * (Xtr_s @ beta)
                sig_neg = 1.0 / (1.0 + np.exp(np.clip(margins, -500, 500)))
                grad    = -(Xtr_s.T @ (y_tr * sig_neg)) / N_tr
                bTz     = float(np.clip(beta @ z_dir, -1e6, 1e6))
                pen     = lam * bTz * z_dir.copy()
                pen[0]  = 0.0
                total   = grad + pen
                # gradient clipping to avoid exploding updates
                gnorm   = np.linalg.norm(total)
                if gnorm > 10.0:
                    total = total * 10.0 / gnorm
                beta   -= lr * total

            logits = Xte_s @ beta
            proba  = 1.0 / (1.0 + np.exp(-logits))
            preds  = np.where(proba >= 0.5, 1.0, -1.0)

            acc_freq_folds[fold_n-1, i]  = float((preds == y_te).mean())
            proxy_freq_folds[fold_n-1, i] = fairness_proxy_freq(beta, z_dir)

        rt_freq_folds[fold_n-1] = time.time() - t0
        print(f"done  t={rt_freq_folds[fold_n-1]:.1f}s")

        # ---- MFVI (ours) — lambda sweep from saved npz ----
        print(f"    MFVI: loading {n_lam} lambdas from npz...", end=" ", flush=True)
        t0 = time.time()
        for i in range(n_lam):
            b = _to_float(fold_data["best_beta"][i])
            S = _to_float(fold_data["best_Sigma"][i])
            if b is None or S is None:
                continue

            inv_S = np.linalg.inv(S)
            proba = vb_logit_pred(X_te, b, S, inv_S)
            preds = np.where(proba >= 0.5, 1.0, -1.0)

            acc_mfvi_folds[fold_n-1, i]   = float((preds == y_te).mean())
            proxy_mfvi_folds[fold_n-1, i] = fairness_proxy_mfvi(b, S, z_te)

            m = compute_all_metrics(y_te, preds, proba, sens_te)
            mfvi_ext_folds[fold_n-1, i, :] = [m[k] for k in metric_keys]

            proba_tr = vb_logit_pred(X_tr, b, S, inv_S)
            preds_tr = np.where(proba_tr >= 0.5, 1.0, -1.0)
            m_tr = compute_all_metrics(y_tr, preds_tr, proba_tr, sens_tr)
            mfvi_ext_tr_folds[fold_n-1, i, :] = [m_tr[k] for k in metric_keys]

        rt_mfvi_folds[fold_n-1] = time.time() - t0
        print(f"done  t={rt_mfvi_folds[fold_n-1]:.1f}s")

    # ------------------------------------------------------------------ #
    # Aggregate
    # ------------------------------------------------------------------ #
    mean_acc_lr    = np.nanmean(acc_lr_folds,    axis=0)
    mean_acc_freq  = np.nanmean(acc_freq_folds,  axis=0)
    std_acc_freq   = np.nanstd(acc_freq_folds,   axis=0)
    mean_acc_mfvi  = np.nanmean(acc_mfvi_folds,  axis=0)
    std_acc_mfvi   = np.nanstd(acc_mfvi_folds,   axis=0)

    mean_proxy_freq = np.nanmean(proxy_freq_folds, axis=0)
    std_proxy_freq  = np.nanstd(proxy_freq_folds,  axis=0)
    mean_proxy_mfvi = np.nanmean(proxy_mfvi_folds, axis=0)
    std_proxy_mfvi  = np.nanstd(proxy_mfvi_folds,  axis=0)

    mean_mfvi_ext = np.nanmean(mfvi_ext_folds, axis=0)   # (n_lam, n_metrics)
    std_mfvi_ext  = np.nanstd(mfvi_ext_folds,  axis=0)
    mean_mfvi_ext_tr = np.nanmean(mfvi_ext_tr_folds, axis=0)
    std_mfvi_ext_tr  = np.nanstd(mfvi_ext_tr_folds,  axis=0)

    mean_rt = np.array([rt_lr_folds.mean(),
                         rt_freq_folds.mean(),
                         rt_mfvi_folds.mean()])
    std_rt  = np.array([rt_lr_folds.std(),
                         rt_freq_folds.std(),
                         rt_mfvi_folds.std()])

    # ------------------------------------------------------------------ #
    # Save
    # ------------------------------------------------------------------ #
    save_path = results_dir / f"comparison_{tag}.npz"
    np.savez(
        save_path,
        vec_lambda      = vec_lambda,
        metric_keys     = np.array(metric_keys),
        mean_acc_lr     = mean_acc_lr,
        mean_acc_freq   = mean_acc_freq,
        std_acc_freq    = std_acc_freq,
        mean_acc_mfvi   = mean_acc_mfvi,
        std_acc_mfvi    = std_acc_mfvi,
        mean_proxy_freq = mean_proxy_freq,
        std_proxy_freq  = std_proxy_freq,
        mean_proxy_mfvi = mean_proxy_mfvi,
        std_proxy_mfvi  = std_proxy_mfvi,
        acc_freq_folds    = acc_freq_folds,
        acc_mfvi_folds    = acc_mfvi_folds,
        proxy_freq_folds  = proxy_freq_folds,
        proxy_mfvi_folds  = proxy_mfvi_folds,
        mfvi_ext_folds    = mfvi_ext_folds,
        mfvi_ext_tr_folds = mfvi_ext_tr_folds,
        mean_mfvi_ext   = mean_mfvi_ext,
        std_mfvi_ext    = std_mfvi_ext,
        mean_mfvi_ext_tr = mean_mfvi_ext_tr,
        std_mfvi_ext_tr  = std_mfvi_ext_tr,
        mean_rt         = mean_rt,
        std_rt          = std_rt,
    )
    print(f"\n  Saved: {save_path}")

    results = dict(
        vec_lambda      = vec_lambda,
        metric_keys     = metric_keys,
        mean_acc_lr     = mean_acc_lr,
        mean_acc_freq   = mean_acc_freq,
        std_acc_freq    = std_acc_freq,
        mean_acc_mfvi   = mean_acc_mfvi,
        std_acc_mfvi    = std_acc_mfvi,
        mean_proxy_freq = mean_proxy_freq,
        std_proxy_freq  = std_proxy_freq,
        mean_proxy_mfvi = mean_proxy_mfvi,
        std_proxy_mfvi  = std_proxy_mfvi,
        acc_freq_folds    = acc_freq_folds,
        acc_mfvi_folds    = acc_mfvi_folds,
        proxy_freq_folds  = proxy_freq_folds,
        proxy_mfvi_folds  = proxy_mfvi_folds,
        mfvi_ext_folds    = mfvi_ext_folds,
        mfvi_ext_tr_folds = mfvi_ext_tr_folds,
        mean_mfvi_ext   = mean_mfvi_ext,
        std_mfvi_ext    = std_mfvi_ext,
        mean_mfvi_ext_tr = mean_mfvi_ext_tr,
        std_mfvi_ext_tr  = std_mfvi_ext_tr,
        mean_rt         = mean_rt,
        std_rt          = std_rt,
        tag             = tag,
        out_dir         = out_dir,
    )

    plot_comparison(results)
    return results


# ---------------------------------------------------------------------------
# Plotting  (called inline — no Agg backend, works in notebooks)
# ---------------------------------------------------------------------------

def plot_comparison(res):
    vec_lambda      = res["vec_lambda"]
    tag             = res["tag"]
    out_dir         = Path(res["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    metric_keys     = res["metric_keys"]

    colors = {"lr": "#4C72B0", "freq": "#DD8452", "mfvi": "#55A868"}
    lw = 2
    fs = 13

    def savefig(fig, name):
        for ext in ("pdf", "png"):
            fig.savefig(out_dir / f"{name}.{ext}", bbox_inches="tight", dpi=150)

    # ------------------------------------------------------------------ #
    # Figure 1: Trade-off accuracy vs fairness proxy — two panels
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"{tag} — Accuracy vs fairness proxy (lambda sweep)", fontsize=fs+1)

    nfp = res["acc_mfvi_folds"].shape[0]

    # Left: accuracy vs lambda  (thin dotted per-fold curves + thick solid mean)
    ax = axes[0]
    ax.axhline(res["mean_acc_lr"].mean(), color=colors["lr"],
               lw=lw, ls="--", label="Standard LR (no fairness)")
    for f in range(nfp):
        ax.semilogx(vec_lambda, res["acc_freq_folds"][f], color=colors["freq"],
                    lw=0.7, ls=":", alpha=0.6)
        ax.semilogx(vec_lambda, res["acc_mfvi_folds"][f], color=colors["mfvi"],
                    lw=0.7, ls=":", alpha=0.6)
    ax.semilogx(vec_lambda, res["mean_acc_freq"], color=colors["freq"],
                lw=lw+0.5, label="Fair LR (freq)  — proxy: $(\\hat{\\beta}^Tz)^2$")
    ax.semilogx(vec_lambda, res["mean_acc_mfvi"], color=colors["mfvi"],
                lw=lw+0.5, label="MFVI (ours)  — proxy: Eq.(29)")
    ax.set_xlabel(r"$\lambda$", fontsize=16)
    ax.set_ylabel("Accuracy (test)", fontsize=fs)
    ax.legend(fontsize=fs-1)
    ax.grid(True, which="both", ls=":")

    # Right: fairness proxy vs lambda  (thin dotted per-fold curves + thick solid mean)
    ax = axes[1]
    for f in range(nfp):
        ax.loglog(vec_lambda, np.maximum(res["proxy_freq_folds"][f], 1e-18),
                  color=colors["freq"], lw=0.7, ls=":", alpha=0.6)
        ax.loglog(vec_lambda, np.maximum(res["proxy_mfvi_folds"][f], 1e-18),
                  color=colors["mfvi"], lw=0.7, ls=":", alpha=0.6)
    ax.loglog(vec_lambda, res["mean_proxy_freq"], color=colors["freq"],
              lw=lw+0.5, label=r"Fair LR (freq): $(\hat{\beta}^Tz)^2$")
    ax.loglog(vec_lambda, res["mean_proxy_mfvi"], color=colors["mfvi"],
              lw=lw+0.5, label=r"MFVI (ours): $(\beta^{*T}z)^2 + \mathrm{Tr}(z^T\Sigma^* z)$")
    ax.set_xlabel(r"$\lambda$", fontsize=16)
    ax.set_ylabel("Fairness proxy (lower = fairer)", fontsize=fs)
    ax.legend(fontsize=fs-1)
    ax.grid(True, which="both", ls=":")

    fig.tight_layout()
    savefig(fig, f"{tag}_comp_fig1_tradeoff")
    plt.show()

    # ------------------------------------------------------------------ #
    # Figure 2: MFVI extended metrics vs lambda
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    fig.suptitle(f"{tag} — MFVI extended metrics vs λ", fontsize=fs+1)
    axes = axes.flatten()

    has_train = "mean_mfvi_ext_tr" in res
    for ai, key in enumerate(metric_keys):
        ax  = axes[ai]
        idx = metric_keys.index(key)
        mu_te  = res["mean_mfvi_ext"][:, idx]
        std_te = res["std_mfvi_ext"][:, idx]
        if has_train:
            mu_tr  = res["mean_mfvi_ext_tr"][:, idx]
            std_tr = res["std_mfvi_ext_tr"][:, idx]
            ax.fill_between(vec_lambda, mu_tr - std_tr, mu_tr + std_tr,
                            color="r", alpha=0.15)
            ax.fill_between(vec_lambda, mu_te - std_te, mu_te + std_te,
                            color="b", alpha=0.15)
            ax.semilogx(vec_lambda, mu_tr, "r", lw=lw, label="Train Set")
            ax.semilogx(vec_lambda, mu_te, "b", lw=lw, label="Test Set")
        else:
            ax.fill_between(vec_lambda, mu_te - std_te, mu_te + std_te,
                            color="b", alpha=0.15)
            ax.semilogx(vec_lambda, mu_te, "b", lw=lw, label="Test Set")
        ax.set_xlabel(r"$\lambda$", fontsize=12)
        ax.set_ylabel(key, fontsize=12)
        ax.grid(True, which="both", ls=":")
        if ai == 0:
            ax.legend(fontsize=10)

    fig.tight_layout()
    savefig(fig, f"{tag}_comp_fig2_mfvi_extended_metrics")
    plt.show()

    # ------------------------------------------------------------------ #
    # Figure 3: Runtime bar chart
    # ------------------------------------------------------------------ #
    method_names = ["Standard LR", "Fair LR (freq)\n[per fold, full sweep]",
                    "MFVI inference\n[per fold, full sweep]"]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(method_names, res["mean_rt"],
                  yerr=res["std_rt"],
                  color=[colors["lr"], colors["freq"], colors["mfvi"]],
                  alpha=0.85, capsize=5)
    ax.set_ylabel("Runtime per fold (s)", fontsize=fs)
    ax.set_title(f"{tag} — Runtime comparison", fontsize=fs)
    ax.grid(axis="y", ls=":", alpha=0.6)
    fig.tight_layout()
    savefig(fig, f"{tag}_comp_fig3_runtime")
    plt.show()

    print(f"  Figures saved to: {out_dir.resolve()}")