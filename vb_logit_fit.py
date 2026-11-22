"""
vb_logit_fit.py
---------------
Coordinate Ascent Variational Inference (CAVI) for *penalised* Bayesian
logistic regression (Fair MFVI).

The model
---------
Likelihood  : p(y_i | x_i, beta) = sigmoid(y_i * x_i^T beta), y_i in {-1, +1}
Prior       : beta | phi ~ N(0, phi^{-1} I),   phi ~ Gamma(a0, rate=b0)
Penalty     : (lambda / 2) * (beta^T z)^2, using H(beta) from Eq. (22)

Setting lambda = 0 recovers the standard Bayesian logistic regression of
Drugowitsch (2013) / Jaakkola & Jordan (2000).

The CAVI updates use D = X.shape[1], including the intercept. After
initializing the precision with a0/b0 and xi0, the updates are:

  a               = a0 + D / 2
  xi_i^2          = max(0, x_i^T Sigma_old x_i + (x_i^T beta_old)^2)
  b               = b0 + (beta_old^T beta_old + Tr(Sigma_old)) / 2
  Sigma_new^{-1}  = (a/b) I + 2 sum_i Omega(xi_i) x_i x_i^T + lambda z z^T
  beta_new        = Sigma_new X^T y / 2

Here Omega(xi) = tanh(xi/2) / (4 xi), with limit 1/8 at zero. Convergence
is checked using the variational lower bound L.

Public API
----------
vb_logit_fit_sigma(X, y, z, lam, a0, b0, xi0)
    -> (beta, Sigma, inv_Sigma, logdet_Sigma, E_a, L, a_n, b_n, xi,
        vec_a, vec_b)
"""

import math
import warnings
import numpy as np
from utils import logdet, nearest_spd


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _lambda_xi(xi: np.ndarray) -> np.ndarray:
    """
    Jaakkola-Jordan auxiliary function:  lam(xi) = tanh(xi/2) / (4*xi).

    The limit as xi -> 0 is 1/8, handled explicitly to avoid 0/0.

    Parameters
    ----------
    xi : ndarray
        Vector of variational parameters (N,).

    Returns
    -------
    ndarray
        lam(xi), same shape as xi.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.tanh(xi / 2.0) / (4.0 * xi)
    # Fix NaN entries that arise at xi == 0
    out = np.where(np.isnan(out), 1.0 / 8.0, out)
    return out


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------


def vb_logit_fit_sigma(
    X, y, z, lam, a0=1e-2, b0=1e-4, xi0=None,
    max_iter=6000,
    spd_retry_log=None,
):
    """Fit penalized Bayesian logistic regression via CAVI.

    The returned E_a is the variational mean of phi, a_n / b_n.
    """
    if not np.isfinite(a0):
        a0 = 1e-2
    if not np.isfinite(b0):
        b0 = 1e-2

    N, D = X.shape

    a_n = a0 + 0.5 * D
    gammaln_an_an = float(math.lgamma(a_n) + a_n)

    t_w = 0.5 * (X * y[:, np.newaxis]).sum(axis=0)

    z = z.ravel()
    penalty_mat = lam * np.outer(z, z)

    if xi0 is None:
        xi0 = np.zeros(N)
    xi = xi0.copy()
    lam_xi = _lambda_xi(xi)

    E_a = a0 / b0

    inv_Sigma = nearest_spd(
        E_a * np.eye(D) + 2.0 * (X.T * lam_xi) @ X + penalty_mat,
        retry_log=spd_retry_log,
    )
    Sigma = np.linalg.inv(inv_Sigma)
    beta = Sigma @ t_w

    b_n = b0 + 0.5 * (beta @ beta + np.trace(Sigma))

    L_last = (
        -N * np.log(2.0)
        + 0.5 * (beta @ inv_Sigma @ beta - logdet(inv_Sigma))
        - (a_n / b_n) * b0
        - a_n * np.log(b_n)
        + gammaln_an_an
    )

    vec_a = [a0, a_n]
    vec_b = [b0, b_n]

    for i in range(max_iter):
        XSX = (X @ Sigma * X).sum(axis=1)
        Xb  = X @ beta
        xi  = np.sqrt(np.maximum(XSX + Xb ** 2, 0.0))
        lam_xi = _lambda_xi(xi)

        b_n = b0 + 0.5 * (beta @ beta + np.trace(Sigma))
        E_a = a_n / b_n

        inv_Sigma = nearest_spd(
            E_a * np.eye(D) + 2.0 * (X.T * lam_xi) @ X + penalty_mat,
            retry_log=spd_retry_log,
        )
        Sigma = np.linalg.inv(inv_Sigma)
        logdet_Sigma = -logdet(inv_Sigma)

        beta = Sigma @ t_w

        L = (
            -np.sum(np.log1p(np.exp(-xi)))
            + np.sum(lam_xi * xi ** 2)
            + 0.5 * (beta @ inv_Sigma @ beta + logdet_Sigma - xi.sum())
            - E_a * b0
            - a_n * np.log(b_n)
            + gammaln_an_an
        )

        vec_a.append(a_n)
        vec_b.append(b_n)

        if (L_last > L) or (abs(L_last - L) < abs(1e-5 * L)):
            break
        L_last = L
    else:
        warnings.warn(
            "vb_logit_fit_sigma: reached maximum number of iterations "
            f"({max_iter}) without convergence.",
            RuntimeWarning,
        )

    L = L - math.lgamma(a0) + a0 * np.log(b0)

    return (
        beta, Sigma, inv_Sigma, logdet_Sigma, E_a, L,
        a_n, b_n, xi, np.array(vec_a), np.array(vec_b),
    )
