"""
utils.py
--------
Utility functions used throughout the Fair MFVI logistic regression pipeline.

Contents
--------
- logdet(A)       : numerically stable log-determinant of a square matrix.
- nearest_spd(A)  : project an arbitrary square matrix onto the nearest
                    Symmetric Positive Definite (SPD) matrix in the
                    Frobenius norm (Higham 1988).

References
----------
Higham, N.J. (1988). "Computing a nearest symmetric positive semidefinite
matrix." Linear Algebra and its Applications, 103, 103-118.
"""

import numpy as np
from scipy.linalg import lu, cholesky, svd
import warnings

# ---------------------------------------------------------------------------
# Log-determinant
# ---------------------------------------------------------------------------

def logdet(A: np.ndarray, use_chol: bool = False) -> float:
    """
    Compute the logarithm of the determinant of a square matrix.

    Uses LU factorisation by default, which is numerically stable even for
    large matrices where ``np.log(np.linalg.det(A))`` would overflow/underflow.
    For symmetric positive-definite matrices pass ``use_chol=True`` to use
    Cholesky factorisation, which is roughly twice as fast.

    Parameters
    ----------
    A : ndarray, shape (n, n)
        Square matrix of float dtype.
    use_chol : bool, optional
        If True, use Cholesky decomposition (A must be SPD). Default: False.

    Returns
    -------
    float
        log |A|. Returns -inf if A is singular.
    """
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square 2-D array.")

    if use_chol:
        # For SPD matrices: log|A| = 2 * sum(log(diag(chol(A))))
        L = cholesky(A, lower=True)
        return 2.0 * np.sum(np.log(np.diag(L)))
    else:
        # General case via LU: A = P @ L @ U
        # log|A| = log|det(P)| + sum(log|diag(U)|)
        # det(P) is +1 or -1, so log|det(P)| = 0.
        import scipy.linalg as sla
        P, L_lu, U = sla.lu(A)
        du = np.diag(U)
        # sign from permutation matrix
        c = np.linalg.det(P) * np.prod(np.sign(du))
        if c == 0:
            return -np.inf
        return np.log(np.abs(c)) + np.sum(np.log(np.abs(du)))


# ---------------------------------------------------------------------------
# Nearest Symmetric Positive Definite matrix
# ---------------------------------------------------------------------------

def nearest_spd(A: np.ndarray, retry_log: list | None = None) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    if A.ndim != 2 or A.shape[1] != n:
        raise ValueError("A must be a square 2-D array.")

    if n == 1:
        if retry_log is not None:
            retry_log.append(0)
        return np.array([[max(float(A[0, 0]), np.finfo(float).eps)]])

    B = (A + A.T) / 2.0
    _, sigma, Vt = svd(B)
    H = Vt.T @ np.diag(sigma) @ Vt
    Ahat = (B + H) / 2.0
    Ahat = (Ahat + Ahat.T) / 2.0

    k = 0
    spacing = np.spacing(np.linalg.norm(Ahat))
    I = np.eye(n)
    max_retries = 300
    while True:
        try:
            cholesky(Ahat, lower=True)
            break
        except np.linalg.LinAlgError:
            k += 1
            min_eig = np.linalg.eigvalsh(Ahat)[0]
            Ahat += (-min_eig * k ** 2 + spacing) * I
            if k >= max_retries:
                Ahat += (abs(min_eig) + 1e-6) * I
                warnings.warn(
                    f"nearest_spd: alcanzado max_retries={max_retries}; "
                    "se aplico diagonal loading de emergencia.",
                    RuntimeWarning,
                )
                break

    if retry_log is not None:
        retry_log.append(k)

    return Ahat
