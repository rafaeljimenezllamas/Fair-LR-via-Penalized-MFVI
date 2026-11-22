"""Fit and evaluate one cross-validation fold for a grid of fairness penalties.

The same seeded log-normal proposals for (a0, b0) are considered at each
penalty. A candidate is accepted when its exact expected training success
rate exceeds the incumbent's by more than 1e-6. The search saves posterior
parameters, predictive summaries, and stopping-rule diagnostics.
"""

import time
import warnings
from pathlib import Path

import numpy as np

from vb_logit_fit import vb_logit_fit_sigma
from vb_logit_pred import vb_logit_pred


def _success_rate(p_y, y_true, n_sim=1000, seed=0):
    N = len(y_true)
    rng = np.random.default_rng(seed)
    U = rng.random((N, n_sim))
    y_pred = 2 * (U < p_y[:, np.newaxis]).astype(float) - 1
    sr_vector = (y_pred == y_true[:, np.newaxis]).mean(axis=0)
    return (float(sr_vector.mean()),
            float(np.percentile(sr_vector, 97.5)),
            float(np.percentile(sr_vector, 2.5)),
            sr_vector)


def _success_rate_exact(p_y, y_true):
    p_correct = np.where(y_true == 1, p_y, 1.0 - p_y)
    return float(p_correct.mean())


def _find_stop(hist_improved, patience, min_samples):
    streak = 0
    for s, improved in enumerate(hist_improved):
        if improved:
            streak = 0
        else:
            streak += 1
        if s + 1 > min_samples and streak >= patience:
            return s + 1
    return len(hist_improved)


def run_fold(
    X, y, z, XT, yT, zT,
    enforce_fairness_on,
    vec_lambda,
    save_path,
    max_samples=500,
    min_samples=15,
    patience=20,
    patience_grid=(1, 5, 10, 15, 20),
    a0_center=1e-2,
    b0_center=1e-4,
    cv=1.0,
    n_sim=1000,
    seed=17,
    mc_seed=12345,
    live_plot=False,
):
    save_path = Path(save_path)
    t_start = time.time()

    z_cavi = z if enforce_fairness_on == "train" else zT
    n_lam = len(vec_lambda)
    patience_grid = list(patience_grid)
    max_patience = max(patience_grid)
    n_pat = len(patience_grid)

    sr_train    = np.zeros(n_lam)
    sr_test     = np.zeros(n_lam)
    sr_train_ex = np.zeros(n_lam)
    sr_test_ex  = np.zeros(n_lam)
    p97_train   = np.zeros(n_lam)
    p2_5_train  = np.zeros(n_lam)
    p97_test    = np.zeros(n_lam)
    p2_5_test   = np.zeros(n_lam)
    dif_train   = np.zeros(n_lam)
    dif_test    = np.zeros(n_lam)
    n_samples_used = np.zeros(n_lam, dtype=int)

    best_beta   = [None] * n_lam
    best_Sigma  = [None] * n_lam
    best_E_a    = [None] * n_lam
    best_L      = [None] * n_lam
    best_vec_a  = [None] * n_lam
    best_vec_b  = [None] * n_lam
    best_a0     = np.full(n_lam, np.nan)
    best_b0     = np.full(n_lam, np.nan)

    # Stopping-rule sensitivity [n_lam, n_pat]
    sens_sr_train = np.zeros((n_lam, n_pat))
    sens_n_used   = np.zeros((n_lam, n_pat), dtype=int)
    sens_time     = np.zeros((n_lam, n_pat))
    sens_a0       = np.full((n_lam, n_pat), np.nan)
    sens_b0       = np.full((n_lam, n_pat), np.nan)

    if live_plot:
        from IPython.display import clear_output
        import matplotlib.pyplot as plt

    sigma_dex = 1.5
    rng = np.random.default_rng(seed)
    log_a = np.log10(a0_center) + sigma_dex * rng.standard_normal(max_samples)
    log_b = np.log10(b0_center) + sigma_dex * rng.standard_normal(max_samples)
    a_samples = 10.0 ** log_a
    b_samples = 10.0 ** log_b

    for i, lam in enumerate(vec_lambda):
        best_sr = -np.inf
        no_improve_streak = 0
        xi_best_lam = np.zeros(len(y))

        hist_sr       = []
        hist_a0       = []
        hist_b0       = []
        hist_improved = []
        hist_time     = []

        for s in range(max_samples):
            a_ini = float(a_samples[s]); b_ini = float(b_samples[s])
            xi_ini = xi_best_lam.copy()

            t_samp = time.time()
            try:
                (beta, Sigma, inv_Sigma, _, E_a, L, a_n, b_n, xi_out,
                 vec_a_hist, vec_b_hist) = vb_logit_fit_sigma(
                    X, y, z_cavi, lam, a_ini, b_ini, xi_ini)
                dt = time.time() - t_samp
            except Exception as exc:
                dt = time.time() - t_samp
                warnings.warn(f"CAVI failed for lam={lam:.3e}, "
                              f"a0={a_ini:.3e}, b0={b_ini:.3e}: {exc}",
                              RuntimeWarning)
                hist_sr.append(-np.inf); hist_a0.append(a_ini); hist_b0.append(b_ini)
                hist_improved.append(False); hist_time.append(dt)
                no_improve_streak += 1
                if s + 1 > min_samples and no_improve_streak >= max_patience:
                    break
                continue

            p_train = vb_logit_pred(X, beta, Sigma, inv_Sigma)
            sr_ex = _success_rate_exact(p_train, y)
            improved = (sr_ex - best_sr > 1e-6)

            hist_sr.append(sr_ex); hist_a0.append(a_ini); hist_b0.append(b_ini)
            hist_improved.append(improved); hist_time.append(dt)

            if improved:
                best_sr = sr_ex
                no_improve_streak = 0
                xi_best_lam = xi_out.copy()

                p_test = vb_logit_pred(XT, beta, Sigma, inv_Sigma)
                sr_tr_mc, p97, p25, _  = _success_rate(p_train, y,  n_sim, seed=mc_seed)
                sr_te_mc, p97t, p25t, _ = _success_rate(p_test,  yT, n_sim, seed=mc_seed)

                sr_train[i]    = sr_tr_mc
                sr_test[i]     = sr_te_mc
                sr_train_ex[i] = sr_ex
                sr_test_ex[i]  = _success_rate_exact(p_test, yT)
                p97_train[i]   = p97
                p2_5_train[i]  = p25
                p97_test[i]    = p97t
                p2_5_test[i]   = p25t
                dif_train[i]   = float((beta @ z) ** 2 + z @ Sigma @ z)
                dif_test[i]    = float((beta @ zT) ** 2 + zT @ Sigma @ zT)
                best_beta[i]   = beta.copy()
                best_Sigma[i]  = Sigma.copy()
                best_E_a[i]    = E_a
                best_L[i]      = L
                best_vec_a[i]  = vec_a_hist.copy()
                best_vec_b[i]  = vec_b_hist.copy()
                best_a0[i]     = a_ini
                best_b0[i]     = b_ini
            else:
                no_improve_streak += 1

            if s + 1 > min_samples and no_improve_streak >= max_patience:
                break

        n_samples_used[i] = len(hist_sr)

        hist_sr_arr   = np.array(hist_sr)
        hist_time_arr = np.array(hist_time)
        for pj, pat in enumerate(patience_grid):
            stop_at = _find_stop(hist_improved, pat, min_samples)
            sub = hist_sr_arr[:stop_at]
            masked = np.where(np.isfinite(sub), sub, -np.inf)
            if np.isfinite(masked).any():
                bidx = int(np.argmax(masked))
                sens_sr_train[i, pj] = sub[bidx]
                sens_a0[i, pj]       = hist_a0[bidx]
                sens_b0[i, pj]       = hist_b0[bidx]
            sens_n_used[i, pj] = stop_at
            sens_time[i, pj]   = hist_time_arr[:stop_at].sum()

        print(f"  lam {i+1:2d}/{n_lam} ({lam:.2e})  |  samples: {len(hist_sr):4d}  "
              f"|  SR exact: {sr_train_ex[i]:.4f}", flush=True)

        if live_plot:
            clear_output(wait=True)
            fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
            m = slice(0, i + 1)
            lam_done = vec_lambda[m]
            axes[0].fill_between(lam_done, p2_5_test[m], p97_test[m], color="b", alpha=0.15)
            axes[0].fill_between(lam_done, p2_5_train[m], p97_train[m], color="r", alpha=0.15)
            axes[0].semilogx(lam_done, sr_train_ex[m], "r", marker="o", ms=3, label="Train")
            axes[0].semilogx(lam_done, sr_test_ex[m],  "b", marker="o", ms=3, label="Test")
            axes[0].set_xlabel(r"$\lambda$"); axes[0].set_ylabel("SR (exact)")
            axes[0].set_title(f"SR vs λ ({i+1}/{n_lam})"); axes[0].legend(); axes[0].grid(True, which="both", ls=":")
            mtr = dif_train[m] > 0; mte = dif_test[m] > 0
            axes[1].loglog(lam_done[mtr], dif_train[m][mtr], "--r", marker="o", ms=3)
            axes[1].loglog(lam_done[mte], dif_test[m][mte],  "--b", marker="o", ms=3)
            axes[1].set_xlabel(r"$\lambda$"); axes[1].set_ylabel(r"$\mathcal{F}$"); axes[1].grid(True, which="both", ls=":")
            Lvals = np.array([bl if bl is not None else np.nan for bl in best_L[:i+1]])
            Fvals = dif_test[m]; valid = (Fvals > 0) & np.isfinite(Lvals)
            if valid.sum() > 1:
                order = np.argsort(Fvals[valid])
                axes[2].semilogx(Fvals[valid][order], Lvals[valid][order], "k", marker="o", ms=3)
            axes[2].set_xlabel(r"$\mathcal{F}$"); axes[2].set_ylabel(r"$\mathcal{L}$"); axes[2].grid(True, which="both", ls=":")
            plt.tight_layout(); plt.show()

    save_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        save_path,
        vec_lambda     = vec_lambda,
        sr_train       = sr_train,
        sr_test        = sr_test,
        sr_train_exact = sr_train_ex,
        sr_test_exact  = sr_test_ex,
        p97_train      = p97_train,
        p2_5_train     = p2_5_train,
        p97_test       = p97_test,
        p2_5_test      = p2_5_test,
        dif_train      = dif_train,
        dif_test       = dif_test,
        n_samples_used = n_samples_used,
        best_a0        = best_a0,
        best_b0        = best_b0,
        patience_grid  = np.array(patience_grid),
        sens_sr_train  = sens_sr_train,
        sens_n_used    = sens_n_used,
        sens_time      = sens_time,
        sens_a0        = sens_a0,
        sens_b0        = sens_b0,
        best_beta      = np.array(best_beta,  dtype=object),
        best_Sigma     = np.array(best_Sigma, dtype=object),
        best_E_a       = np.array(best_E_a),
        best_L         = np.array(best_L),
        best_vec_a     = np.array(best_vec_a, dtype=object),
        best_vec_b     = np.array(best_vec_b, dtype=object),
        y = y, yT = yT, z = z, zT = zT,
    )

    elapsed = time.time() - t_start
    print(f"\n  Fold saved to {save_path}  (total: {elapsed/60:.1f} min)")
