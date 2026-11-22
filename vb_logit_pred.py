"""
vb_logit_pred.py
----------------
Posterior predictive probabilities for a fitted Bayesian logistic regression.

For each new observation x (a row of X), the function approximates

    p(y = 1 | x, X_train, y_train)

using the variational posterior N(beta | w, V) returned by
``vb_logit_fit_sigma``.  The computation follows Eqs. (13)-(16) of the
paper and the original derivation in Drugowitsch (2013) / Jaakkola &
Jordan (2000).

Public API
----------
vb_logit_pred(X, w, V, inv_V) -> ndarray, shape (K,)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _lambda_xi(xi: np.ndarray) -> np.ndarray:
    """
    Jaakkola-Jordan function  lam(xi) = tanh(xi/2) / (4*xi),
    with the limit 1/8 at xi = 0.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.tanh(xi / 2.0) / (4.0 * xi)
    return np.where(np.isnan(out), 1.0 / 8.0, out)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def vb_logit_pred(
    X: np.ndarray,
    w: np.ndarray,
    V: np.ndarray,
    inv_V: np.ndarray,
    max_iter: int = 1000,
) -> np.ndarray:
    """
    Compute p(y = 1 | x, X_train, y_train) for each row x of X.

    The posterior N(beta | w, V) is used to integrate out beta analytically,
    following Eqs. (14)-(16) and the closed-form result (13) of the paper.

    Each test point is handled independently, iterating on the variational
    parameter xi_i until the per-point lower bound converges.

    Parameters
    ----------
    X : ndarray, shape (K, D)
        K test observations, one per row.
    w : ndarray, shape (D,)
        Posterior weight mean (beta*).
    V : ndarray, shape (D, D)
        Posterior weight covariance (Sigma*).
    inv_V : ndarray, shape (D, D)
        Inverse of V (Sigma*^{-1}).  Passed explicitly to avoid recomputing.
    max_iter : int, optional
        Maximum number of xi update iterations per test point. Default: 1000.

    Returns
    -------
    ndarray, shape (K,)
        Predictive probability p(y = 1 | x_k) for each test point k.
    """
    N, Dx = X.shape

    # ------------------------------------------------------------------ #
    # Pre-compute quantities that are re-used across iterations
    # (vectorised over all K test points simultaneously)
    #
    # Notation follows Drugowitsch (2013) and Eqs. (14)-(16):
    #   w_t   = V^{-1} w + x/2          (D-vector per test point)
    #   Vx    = V x                      (D-vector per test point)
    #   xVx   = x^T V x                  (scalar per test point)
    # ------------------------------------------------------------------ #

    # w_t[i] = inv_V @ w + 0.5 * x_i    ->  shape (K, D)
    w_t      = np.ones((N, 1)) * (w @ inv_V) + 0.5 * X   # broadcast
    Vx       = X @ V                                       # (K, D)
    xVx      = (Vx * X).sum(axis=1)                       # (K,)  x^T V x
    xVx2     = xVx ** 2                                    # (K,)

    # V x x^T V^T w_t  ->  shape (K, D)
    #   = Vx * ( (w_t * Vx).sum(axis=1) )[:, None]
    VxxVwt   = Vx * (np.einsum("id,id->i", w_t, Vx))[:, np.newaxis]

    # V w_t  ->  shape (K, D)
    Vwt      = w_t @ V

    # ------------------------------------------------------------------ #
    # Initialisation: xi = 0  =>  lam_xi = 1/8
    # (Sherman-Morrison / matrix-determinant lemma for xi -> 0 case)
    # ------------------------------------------------------------------ #
    xi       = np.zeros(N)
    lam_xi   = _lambda_xi(xi)                              # = 1/8 everywhere

    # a_xi = 1 / (4 + x^T V x)
    a_xi     = 1.0 / (4.0 + xVx)                          # (K,)

    # w_xi = V w_t - a_xi * V x x^T V^T w_t
    w_xi     = Vwt - a_xi[:, np.newaxis] * VxxVwt         # (K, D)

    # log|V_xi| = -log(1 + x^T V x / 4)
    logdetV_xi = -np.log1p(xVx / 4.0)                     # (K,)

    # w_xi^T inv_V w_xi + (w_xi^T x)^2 / 4
    wVw_xi   = (w_xi @ inv_V * w_xi).sum(axis=1) + (w_xi * X).sum(axis=1) ** 2 / 4.0

    L_last   = 0.5 * (logdetV_xi.sum() + wVw_xi.sum()) - N * np.log(2.0)

    # ------------------------------------------------------------------ #
    # Iterate the xi updates (EM on the per-point variational bound)
    # ------------------------------------------------------------------ #
    for _ in range(max_iter):
        # -- Eq. 16: update xi --
        xi     = np.sqrt(np.maximum(
            xVx - a_xi * xVx2 + (w_xi * X).sum(axis=1) ** 2, 0.0
        ))
        lam_xi = _lambda_xi(xi)

        # -- Eqs. 14-15: Sherman-Morrison update of Sigma_hat and beta_hat --
        a_xi       = 2.0 * lam_xi / (1.0 + 2.0 * lam_xi * xVx)   # (K,)
        w_xi       = Vwt - a_xi[:, np.newaxis] * VxxVwt            # (K, D)
        logdetV_xi = -np.log1p(2.0 * lam_xi * xVx)                # (K,)

        # Per-point quadratic form in the bound
        wVw_xi = (
            (w_xi @ inv_V * w_xi).sum(axis=1)
            + 2.0 * lam_xi * (w_xi * X).sum(axis=1) ** 2
        )

        # -- Variational lower bound (Eq. 13), constant terms omitted --
        L = np.sum(
            0.5 * (logdetV_xi + wVw_xi - xi)
            - np.log1p(np.exp(-xi))
            + lam_xi * xi ** 2
        )

        # Bound must not decrease (stop if it does – numerical noise)
        if L_last > L:
            break
        if abs(L_last - L) < abs(1e-5 * L):
            break
        L_last = L

    # ------------------------------------------------------------------ #
    # Closed-form predictive probability (Eq. 13)
    # ------------------------------------------------------------------ #
    # p(y=1 | x) = sigma(xi) / sqrt(1 + 2*lam_xi*xVx)
    #              * exp( 0.5*(-xi - w^T inv_V w + w_xi_Vw_xi) + lam_xi*xi^2 )
    w_inv_V_w = float(w @ inv_V @ w)   # scalar, shared across all test points

    log_pred = (
        -np.log1p(np.exp(-xi))          # log sigma(xi)
        - 0.5 * np.log1p(2.0 * lam_xi * xVx)
        + 0.5 * (-xi - w_inv_V_w + wVw_xi)
        + lam_xi * xi ** 2
    )
    out = np.exp(log_pred)

    # Clip to [0, 1] for safety (should already be in range)
    return np.clip(out, 0.0, 1.0)
