"""Numerical utilities used by the scVOLT encoder and analysis helpers."""

from typing import List
import torch
from scipy.spatial.distance import cdist
import numpy as np
from scipy import stats

def fdr_adjust(pvals, method="BH"):
    """Adjust p-values with Benjamini-Hochberg or Benjamini-Yekutieli FDR control.

    Args:
        pvals (array-like): Raw p-values. Non-finite values are preserved as NaN.
        method (str): FDR method, either ``"BH"`` or ``"BY"``.

    Returns:
        np.ndarray: Adjusted q-values in the original input order.
    """
    pvals = np.asarray(pvals, dtype=float)
    qvals = np.full_like(pvals, np.nan, dtype=float)

    # Keep only finite values for ranking while preserving missing entries.
    ok = np.isfinite(pvals)
    pv = pvals[ok]
    m = pv.size
    if m == 0:
        return qvals

    order = np.argsort(pv)
    pv_sorted = pv[order]

    method = method.upper()
    if method == "BH":
        scale = 1.0
    elif method == "BY":
        scale = np.sum(1.0 / np.arange(1, m + 1))
    else:
        raise ValueError("method must be 'BH' or 'BY'")

    # Compute ranked q-values and enforce monotonicity from largest to smallest p-value.
    denom = np.arange(1, m + 1)
    q_sorted = pv_sorted * (m * scale) / denom
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    out = np.empty_like(q_sorted)
    out[order] = q_sorted
    qvals[ok] = out
    return qvals

def pearson_pvals_vectorized(X, y):
    """Compute Pearson correlations and two-sided p-values for many features.

    Args:
        X (array-like): Feature matrix with shape ``(n_samples, n_features)``.
        y (array-like): Response vector with shape ``(n_samples,)``.

    Returns:
        tuple[np.ndarray, np.ndarray]: Correlation coefficients and p-values.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    n, m = X.shape
    if y.shape[0] != n:
        raise ValueError("X and y must have the same number of samples")

    # Center the response and each feature before computing covariance.
    y0 = y - np.nanmean(y)
    X0 = X - np.nanmean(X, axis=0, keepdims=True)

    # Compute sample standard deviations with NaN-tolerant reductions.
    y_std = np.sqrt(np.nansum(y0 * y0) / max(n - 1, 1))
    X_std = np.sqrt(np.nansum(X0 * X0, axis=0) / max(n - 1, 1))

    r = np.zeros(m, dtype=float)
    p = np.ones(m, dtype=float)

    # Constant responses or constant feature columns keep the default r=0 and p=1.
    valid = (y_std > 0) & (X_std > 0) & np.isfinite(X_std)
    if not valid.any():
        return r, p

    # Convert covariance to Pearson r and then to the standard t statistic.
    cov = np.nansum(X0[:, valid] * y0[:, None], axis=0) / max(n - 1, 1)
    rv = cov / (X_std[valid] * y_std)
    rv = np.clip(rv, -0.999999999, 0.999999999)

    df = n - 2
    tv = rv * np.sqrt(df / (1.0 - rv * rv))
    pv = 2.0 * stats.t.sf(np.abs(tv), df=df)

    r[valid] = rv
    p[valid] = pv
    return r, p


def reference_dataset(
    X, dtype: torch.dtype, device: torch.device
) -> torch.Tensor:
    """Select features, transpose dataset and convert to Tensor.

    Args:
        X (array-like): The input data
        dtype (torch.dtype): The dtype to create
        device (torch.device): The device to create on
        keep_idx (Iterable): The variables to keep.

    Returns:
        torch.Tensor: The reference dataset A.
    """

    # Transpose the matrix so features are rows and samples are columns.
    A = X.T

    # Check that the dataset is positive.
    if not (A >= 0).all():
        raise ValueError("Input expression values must be non-negative")

    # If the dataset is sparse, make it dense.
    try:
        A = A.todense()
    except AttributeError:
        pass

    # Send the matrix `A` to PyTorch.
    return torch.from_numpy(A).to(device=device, dtype=dtype).contiguous()


def compute_ground_cost(
    features,
    cost: str,
    eps: float,
    dtype: torch.dtype,
    device: torch.device,
) -> torch.Tensor:
    """Compute the dense ground-cost kernel.

    Args:
        features (array-like): An array with the features to compute the cost on.
        cost (str): The function to compute the cost. Scipy distances are allowed.
        eps (float): Entropic regularization used to scale the distances.
        dtype (torch.dtype): The dtype for the output.
        device (torch.device): The device for the output.

    Returns:
        torch.Tensor: The ground-cost kernel.
    """

    if eps <= 0:
        raise ValueError("eps must be positive")

    # Compute pairwise distances and scale them into an entropic kernel.
    K = cdist(features, features, metric=cost)
    if not np.isfinite(K).all():
        raise ValueError("Ground-cost matrix contains non-finite values")
    max_cost = K.max()
    if max_cost == 0:
        return torch.ones(K.shape, device=device, dtype=dtype)

    K = torch.from_numpy(K).to(device=device, dtype=dtype)
    K /= eps * max_cost
    K = torch.exp(-K).to(device=device, dtype=dtype)
    return K


def normalize_tensor(X: torch.Tensor) -> torch.Tensor:
    """Normalize a tensor along columns.

    Args:
        X (torch.Tensor): The tensor to normalize.

    Returns:
        torch.Tensor: The normalized tensor.
    """

    return X / X.sum(0)


def entropy(
    X: torch.Tensor, min_one: bool = False, rescale: bool = False
) -> torch.Tensor:
    """Entropy function, :math:`E(X) = \langle X, \log X - 1 \rangle`.

    Args:
        X (torch.Tensor):
            The parameter to compute the entropy of.
        min_one (bool, optional):
            Whether to inclue the :math:`-1` in the formula. Defaults to False.
        rescale (bool, optional):
            Rescale so that the value is between 0 and 1 (when min_one=False).

    Returns:
        torch.Tensor: The entropy of X.
    """
    offset = 1 if min_one else 0
    scale = X.shape[1] * np.log(X.shape[0]) if rescale else 1
    return -torch.sum(X * (torch.nan_to_num(X.log()) - offset)) / scale


def entropy_dual_loss(Y: torch.Tensor) -> torch.Tensor:
    """Compute the Legendre dual of the entropy.

    Args:
        Y (torch.Tensor): The input parameter.

    Returns:
        torch.Tensor: The loss.
    """
    return -torch.logsumexp(Y, dim=0).sum()


def ot_dual_loss(
    A: dict, G: dict, K: dict, eps: float, gene_weights: torch.Tensor, dim=(0, 1)
) -> torch.Tensor:
    """Compute the Legendre dual of the entropic OT loss.

    Args:
        A (dict): The input data.
        G (dict): The dual variable.
        K (dict): The kernel.
        eps (float): The entropic regularization.
        gene_weights (torch.Tensor): The per-gene weights.
        dim (tuple, optional): How to sum the loss. Defaults to (0, 1).

    Returns:
        torch.Tensor: The loss
    """
    log_fG = G / eps

    # Stabilize the matrix product by subtracting the column-wise maximum first.
    scale = log_fG.max(dim=0).values
    prod = torch.log(K @ torch.exp(log_fG - scale)) + scale
    # Weight the stabilized product by the input distribution and gene weights.
    return eps * torch.sum(gene_weights * A * prod, dim=dim)


def early_stop(history: List, tol: float, nonincreasing: bool = False) -> bool:
    """Based on a history and a tolerance, whether to stop early or not.

    Args:
        history (List):
            The loss history.
        tol (float):
            The tolerance before early stopping.
        nonincreasing (bool, optional):
            When False, throws an error if the loss goes up. Defaults to False.

    Raises:
        ValueError: When the loss goes up.

    Returns:
        bool: Whether to stop early.
    """
    # Stop immediately on invalid losses so the caller can surface the failure.
    if len(history) > 0 and not torch.isfinite(history[-1]):
        raise ValueError("Error: Loss is not finite!")

    # Require a short warm-up before comparing adjacent losses.
    if len(history) < 3:
        return False

    # Optionally stop when the monitored loss increases beyond the tolerance.
    if nonincreasing and (history[-1] - history[-3]) > tol:
        return True

    # Stop once successive losses are within tolerance.
    if abs(history[-1] - history[-2]) < tol:
        return True

    # Continue optimization otherwise.
    return False
