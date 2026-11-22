"""Draw seeded posterior histograms of H(beta) on logarithmic axes."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter


def _sample_H(beta, Sigma, z, n_sim, rng, floor):
    S = (Sigma + Sigma.T) / 2.0            # Symmetrize before sampling.
    samples = rng.multivariate_normal(beta, S, size=n_sim)
    return np.clip((samples @ z) ** 2, floor, None)


def _plot_hist_logx(ax, H, color, label, edges):
    """Draw the log10(H) density as bars on a logarithmic x-axis."""
    counts, _ = np.histogram(np.log10(H), bins=edges, density=True)
    lefts  = 10.0 ** edges[:-1]
    widths = 10.0 ** edges[1:] - 10.0 ** edges[:-1]
    ax.bar(lefts, counts, width=widths, align="edge",
           alpha=0.55, color=color, label=label, edgecolor="none")


def _set_log_xticks(ax, edges):
    """Choose log-axis ticks for the displayed range."""
    xmin = 10.0 ** edges[0]
    xmax = 10.0 ** edges[-1]
    R = edges[-1] - edges[0]   # Width in decades.

    # Use one major label per decade for wide ranges.
    if R >= 2.5:
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))

        def _fmt_decade(x, pos=None):
            if x <= 0:
                return ""
            k = int(round(np.log10(x)))
            if np.isclose(x, 10.0 ** k):
                return rf"$10^{{{k}}}$"
            return ""

        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_decade))
        ax.xaxis.set_minor_locator(
            LogLocator(base=10, subs=tuple(np.arange(2, 10) * 0.1))
        )
        ax.xaxis.set_minor_formatter(NullFormatter())
        return

    # Use at most four roughly log-spaced labels for narrow ranges.
    nice = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
    kmin = int(np.floor(np.log10(xmin))) - 1
    kmax = int(np.ceil(np.log10(xmax))) + 1
    cand = []
    for k in range(kmin, kmax + 1):
        decade = 10.0 ** k
        for m in nice:
            x = m * decade
            if xmin <= x <= xmax:
                cand.append(x)
    cand = np.array(sorted(set(cand)), dtype=float)

    if len(cand) > 4:
        targets = np.linspace(np.log10(xmin), np.log10(xmax), 4)
        chosen = []
        for t in targets:
            j = int(np.argmin(np.abs(np.log10(cand) - t)))
            chosen.append(cand[j])
        ticks = np.array(sorted(set(chosen)), dtype=float)

        # Fill duplicate selections with the most distant remaining candidates.
        while len(ticks) < min(4, len(cand)):
            remaining = np.array([x for x in cand if x not in ticks])
            if len(remaining) == 0:
                break
            d = np.array([
                np.min(np.abs(np.log10(x) - np.log10(ticks)))
                for x in remaining
            ])
            ticks = np.sort(np.append(ticks, remaining[np.argmax(d)]))
    else:
        ticks = cand

    ax.set_xticks(ticks)

    def _fmt_narrow(x, pos=None):
        # Use decimals for values between 0.001 and 1.
        if 1e-3 <= x < 1.0:
            return f"{x:.3g}"
        exponent = int(np.floor(np.log10(x)))
        mantissa = x / (10.0 ** exponent)
        if np.isclose(mantissa, 1.0):
            return rf"$10^{{{exponent}}}$"
        return rf"${mantissa:g}\times10^{{{exponent}}}$"

    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_narrow))
    ax.xaxis.set_minor_locator(
        LogLocator(base=10, subs=tuple(np.arange(2, 10) * 0.1))
    )
    ax.xaxis.set_minor_formatter(NullFormatter())


def h_beta_histograms(
    tag,
    lam_indices=None,          # Paper grid indices (0, 32, 59).
    lam_values=None,           # Alternatively, select the nearest values.
    fold=1,
    which="test",              # "test", "train", or "both".
    n_sim=20000,
    n_bins=40,
    floor=1e-13,
    seed=0,
    layout="col",              # Arrange panels in a column or a row.
    share_bins=False,          # Use a separate x-axis range for each panel by default.
    results_dir="results",
    figures_dir="figures",
    show=True,
):
    """
    Draw one panel per requested penalty using lam_indices or lam_values.

    With share_bins=True, all panels use the same x-axis limits.
    """
    results_dir = Path(results_dir)
    out_dir = Path(figures_dir) / tag.lower()
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    data = np.load(results_dir / f"Results_{tag}_CV_{fold}.npz", allow_pickle=True)
    vec_lambda = np.asarray(data["vec_lambda"], dtype=float)
    betas  = data["best_beta"]
    sigmas = data["best_Sigma"]
    z  = np.asarray(data["z"],  dtype=float).ravel()
    zT = np.asarray(data["zT"], dtype=float).ravel()

    # Select penalty indices.
    if lam_indices is None and lam_values is None:
        lam_indices = (0, len(vec_lambda) // 2, len(vec_lambda) - 1)
    if lam_values is not None:
        lam_indices = tuple(int(np.argmin(np.abs(vec_lambda - v))) for v in lam_values)

    # Sample H for each panel.
    H_by_panel = {}
    for idx in lam_indices:
        b, S = betas[idx], sigmas[idx]
        if b is None or S is None:
            H_by_panel[idx] = None
            continue
        b = np.asarray(b, dtype=float); S = np.asarray(S, dtype=float)
        entry = {}
        if which in ("train", "both"):
            entry["train"] = _sample_H(b, S, z,  n_sim, rng, floor)
        if which in ("test", "both"):
            entry["test"]  = _sample_H(b, S, zT, n_sim, rng, floor)
        H_by_panel[idx] = entry

    def _edges_from(logvals):
        lo, hi = logvals.min(), logvals.max()
        pad = 0.05 * (hi - lo) if hi > lo else 0.5
        return np.linspace(lo - pad, hi + pad, n_bins)

    shared_edges = None
    if share_bins:
        all_logH = np.concatenate([
            np.log10(arr) for e in H_by_panel.values() if e for arr in e.values()
        ])
        shared_edges = _edges_from(all_logH)

    # ---- 2) plot ----
    idxs = list(lam_indices)
    if layout == "col":
        fig, axes = plt.subplots(len(idxs), 1, figsize=(7, 3.4 * len(idxs)),
                                 squeeze=False)
        axes = axes[:, 0]
    else:
        fig, axes = plt.subplots(1, len(idxs), figsize=(6 * len(idxs), 4.2),
                                 squeeze=False)
        axes = axes[0]

    for ax, idx in zip(axes, idxs):
        entry = H_by_panel[idx]
        lam_val = vec_lambda[idx]
        if entry is None:
            ax.set_title(fr"$\lambda={lam_val:.1e}$  (no data)")
            continue

        if shared_edges is not None:
            edges = shared_edges
        else:
            logvals = np.concatenate([np.log10(a) for a in entry.values()])
            edges = _edges_from(logvals)

        if "train" in entry:
            _plot_hist_logx(ax, entry["train"], "#C44E52", "Train Set", edges)
        if "test" in entry:
            _plot_hist_logx(ax, entry["test"], "#6b6be0", "Test Set", edges)

        ax.set_xscale("log")
        ax.set_xlim(10.0 ** edges[0], 10.0 ** edges[-1])
        _set_log_xticks(ax, edges)
        ax.set_xlabel(r"$H(\beta)$", fontsize=18)
        ax.set_ylabel(r"$f(\cdot)$", fontsize=18)
        ax.set_title(fr"$\lambda = {lam_val:.1e}$", fontsize=13)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=12)
        ax.grid(True, which="major", ls=":", alpha=0.45)
        ax.grid(True, which="minor", axis="x", ls=":", alpha=0.15)

    fig.tight_layout()
    itag = "_".join(str(i) for i in idxs)
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{tag}_H_beta_hist_{itag}.{ext}",
                    dpi=150, bbox_inches="tight")
    if show:
        plt.show()
        plt.close(fig)
        return None
    return fig