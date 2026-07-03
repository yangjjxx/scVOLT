"""Analysis and reproducibility helpers for scVOLT workflows."""
from typing import Optional
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
import random
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
import torch


def quality_control(input_data: dict) -> None:
    """Filter low-quality cells and weakly detected genes in-place.

    Args:
        input_data (dict): Mapping from dataset key to expression matrix.
    """

    for key in input_data.keys():
        # Calculate mitochondrial gene ratio for each cell.
        mt_genes = input_data[key].columns.str.startswith("MT-")
        mt_genes_ratio = input_data[key].apply(lambda row: sum(row[mt_genes]) / sum(row), axis=1)
        # Count how many cells express each gene.
        gene_detection_counts = input_data[key].apply(lambda col: sum(1 for i in col if i != 0), axis=0)
        # Count how many genes are detected in each cell.
        genes_per_cell = input_data[key].apply(lambda row: sum(1 for i in row if i != 0), axis=1)
        # Drop low-quality cells and rarely detected genes.
        filtered_cells = input_data[key].loc[(genes_per_cell <= 200) | (mt_genes_ratio >= 0.3), :].index
        filtered_genes = input_data[key].loc[:, gene_detection_counts <= 3].columns
        input_data[key] = input_data[key].drop(labels=filtered_cells, axis=0)
        input_data[key] = input_data[key].drop(labels=filtered_genes, axis=1)


def get_metrics(Y_true, Y_pred) -> dict:
    """Compute binary classification metrics and effect size.

    Args:
        Y_true (array-like): Binary ground-truth labels.
        Y_pred (array-like): Predicted scores or probabilities.

    Returns:
        dict: AUROC, AUPRC, and Cohen's d for the positive versus negative groups.
    """
    Y_true = np.asarray(Y_true).astype(int).reshape(-1)
    Y_pred = np.asarray(Y_pred).astype(float).reshape(-1)

    # Compute threshold-free classification metrics from continuous scores.
    auroc = roc_auc_score(Y_true, Y_pred)
    auprc = average_precision_score(Y_true, Y_pred)

    # Split scores by class before estimating Cohen's d.
    group1 = Y_pred[Y_true == 1]
    group0 = Y_pred[Y_true == 0]

    cohens_d = np.nan

    if group1.size >= 2 and group0.size >= 2:
        n1, n0 = group1.size, group0.size
        s1 = np.std(group1, ddof=1)
        s0 = np.std(group0, ddof=1)

        pooled_std = np.sqrt(
            ((n1 - 1) * s1**2 + (n0 - 1) * s0**2) / (n1 + n0 - 2)
        )

        if pooled_std > 0:
            cohens_d = float((np.mean(group1) - np.mean(group0)) / pooled_std)

    return {
        "AUROC": float(auroc),
        "AUPRC": float(auprc),
        "Cohens_d": cohens_d,
    }


def set_model_seed(random_seed: Optional[int] = 100) -> None:
    """Set random seeds for NumPy, PyTorch, and Python's random module.

    Args:
        random_seed (int|None): Seed value. If ``None``, no seed is set.
    """
    if random_seed is not None:
        # Seed CPU-side libraries first.
        np.random.seed(random_seed)
        torch.manual_seed(random_seed)
        random.seed(random_seed)

        if torch.cuda.is_available():
            # Seed all CUDA devices and request deterministic kernels where available.
            torch.cuda.manual_seed(random_seed)
            torch.cuda.manual_seed_all(random_seed)  
            torch.use_deterministic_algorithms(True)

