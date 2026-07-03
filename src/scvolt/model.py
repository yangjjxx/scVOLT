"""High-level scVOLT model API for encoding data and predicting drug sensitivity."""

import torch
import numpy as np
import pandas as pd

from pyreadr import read_r
from .encoder import MatrixFactorizationEncoder
from .predictor import ElasticNetPredictor
from importlib.resources import files
from scipy.stats import t as t_dist
from scipy.stats import rankdata

pkg_path = files("scvolt")

_DEFAULT_H_REGULARIZATION = {
    "bk": 1e-3,
    "sc": 1e-3,
    "st": 1e-3,
}


class SCVOLTModel():
    """Coordinate scVOLT encoding and ElasticNet-based sensitivity prediction."""

    def __init__(
            self,
            drug_name: str,
            latent_dim: int = 50,
            num_DRGs: int = 1200,
            h_regularization=None,
            w_regularization: float = 1e-3,
            eps: float = 5e-3,
            cost: str = "cosine",
            
    ) -> None:
        """Initialize model hyperparameters and the matrix factorization encoder.

        Args:
            drug_name (str): Drug name matching a CTRP response column and signature file.
            latent_dim (int): Number of latent factors in the shared embedding.
            num_DRGs (int): Number of signature genes to keep.
            h_regularization (dict | float | None): Entropy regularization for modality-specific factors.
            w_regularization (float): Entropy regularization for the shared factor.
            eps (float): Entropic regularization for optimal transport.
            cost (str): SciPy distance metric used to build the transport ground cost.
        """

        if h_regularization is None:
            h_regularization = dict(_DEFAULT_H_REGULARIZATION)

        # Validate entropy regularization values before constructing the encoder.
        if isinstance(h_regularization, dict):
            for resolution in h_regularization:
                if h_regularization[resolution] <= 0:
                    raise ValueError("h_regularization values must be positive")
        else:
            if h_regularization <= 0:
                raise ValueError("h_regularization must be positive")

        self.drug_name = drug_name
        self.num_DRGs = num_DRGs

        self.encoder = MatrixFactorizationEncoder(
            latent_dim=latent_dim,
            w_regularization=w_regularization,
            h_regularization=h_regularization,
            eps=eps,
            cost=cost
        )

    def encode(
            self,
            data: dict,
            max_iter_inner: int = 800,
            max_iter: int = 30,
            device=None,
            dtype: torch.dtype = torch.double,
            lr: float = 1,
            optim_name: str = "lbfgs",
            tol_inner: float = 1e-5,
            tol_outer: float = 1e-3,
            normalize_rows: bool = False,
            enable_gene_weight: bool = False,
            gene_weight: dict = None,
    ) -> dict:
        """Encode one input expression matrix together with the bulk reference data.

        Args:
            data (dict): Single-entry mapping from input key to expression matrix.
            max_iter_inner (int): Maximum optimizer steps for each inner update.
            max_iter (int): Maximum number of alternating encoder updates.
            device (torch.device | str | None): Device used by PyTorch tensors. If None, CUDA is used when available.
            dtype (torch.dtype): Tensor dtype used during optimization.
            lr (float): Optimizer learning rate.
            optim_name (str): Optimizer name, one of ``"lbfgs"``, ``"sgd"``, or ``"adam"``.
            tol_inner (float): Early-stopping tolerance for inner updates.
            tol_outer (float): Early-stopping tolerance for alternating updates.
            normalize_rows (bool): Reserved normalization flag passed to the encoder.
            enable_gene_weight (bool): Whether to pass per-gene weights into the encoder.
            gene_weight (dict): Optional per-resolution gene weights.

        Returns:
            dict: Encoded bulk and input latent factors.
        """
        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            device = torch.device(device)

        self.key, input_data = next(iter(data.items()))
        
        # Load bulk expression and drug-response resources bundled with the package.
        bk_exp = read_r(pkg_path / "resources/Drug Screens/CCLE_Expr_Log2TPM.rds")[None]

        bk_response = read_r(pkg_path / "resources/Drug Screens/CTRP2_Res_nAUC.rds")[None]

        # Align bulk expression and response matrices to their shared samples.
        common_sample = bk_response.index.intersection(bk_exp.index)
        bk_response = bk_response.loc[common_sample, :]
        bk_exp = bk_exp.loc[common_sample, :]

        # Remove genes with missing expression before selecting shared signature genes.
        bk_exp = bk_exp.dropna(axis=1)
        input_data = input_data.dropna(axis=1)

        # Select the requested drug response and its ranked signature genes.
        self.y = bk_response[self.drug_name].dropna()
        sig_genes = pd.read_table(pkg_path / "resources/drgs/ctrp/" / f"{self.drug_name}_signature.txt",
                                  sep=' ')['x'].to_list()
        
        # Keep only signature genes that exist in both bulk and input expression matrices.
        bk_input_genes = set(bk_exp.columns) & set(input_data.columns)
        common_signature = [g for g in sig_genes if g in bk_input_genes]

        if self.num_DRGs > len(common_signature):
            raise ValueError(
                f"Requested {self.num_DRGs} genes, but only "
                f"{len(common_signature)} signature genes are available"
            )
        common_vgenes = common_signature[0:self.num_DRGs]

        input_data = input_data.loc[:, common_vgenes]

        # Build the aligned bulk matrix for the same selected signature genes.
        bk_data = bk_exp.loc[:, common_vgenes]
        bk_data = bk_data.reindex(self.y.index)

        rna_data={'bk': bk_data, self.key: input_data}
        
        # Train the encoder and return all latent factors.
        encoded_data = self.encoder.train(
            rna_data=rna_data,
            max_iter_inner=max_iter_inner,
            max_iter=max_iter,
            device=device,
            dtype=dtype,
            lr=lr,
            optim_name=optim_name,
            tol_inner=tol_inner,
            tol_outer=tol_outer,
            normalize_rows=normalize_rows,
            enable_gene_weight=enable_gene_weight,
            gene_weight=gene_weight
        )
        
        return encoded_data

    def predict(
            self,
            encoded_data: dict, 
            max_iter: int = 5000,
            control_fdr: bool = False,
            alpha: float = 0.05
            ) -> dict:
        """Predict drug sensitivity scores from encoded latent factors.

        Args:
            encoded_data (dict): Encoder output containing ``"H_bk"`` and the input key.
            max_iter (int): Maximum solver iterations for ElasticNetCV.
            control_fdr (bool): Whether to filter latent factors by Spearman FDR.
            alpha (float): FDR threshold used when ``control_fdr`` is enabled.

        Returns:
            dict: Bulk and input sensitivity scores.
        """

        predictor = ElasticNetPredictor(max_iter=max_iter)
        encoded_bk  = encoded_data["H_bk"]
        encoded_input = encoded_data["H_"+self.key]

        if control_fdr:

            # Select latent factors associated with the bulk drug response.
            pvals = spearman_r_pvals_vectorized(encoded_bk, self.y)
            qvals = fdr_adjust(pvals, method="BH")
            idx = np.where(qvals <= alpha)[0]

            if idx.size == 0:
                print("No features passed FDR threshold; using all dims.")
            else:
                # Keep the same selected latent dimensions for training and prediction.
                encoded_bk = encoded_bk[:, idx]
                encoded_input = encoded_input[:, idx]
                print("selected dims:", encoded_bk.shape[1])

            predictor.train(encoded_bk, self.y)
            bk_pre = predictor.predict(encoded_bk)
            input_pre = predictor.predict(encoded_input)
            results = {"bk_scores": bk_pre,
                       self.key+"_scores": input_pre}
         
        else:
            # Train on bulk factors and apply the same predictor to the input factors.
            predictor.train(encoded_bk, self.y)
            bk_pre = predictor.predict(encoded_bk)
            input_pre = predictor.predict(encoded_input)

            results = {"bk_scores": bk_pre,
                       self.key+"_scores": input_pre}
        
        return results
                    

    @staticmethod
    def get_drugs() -> list:
        """Return all available drug names from the bundled CTRP response table.

        Returns:
            list: Drug names available for scVOLT prediction.
        """
        bk_response = read_r(pkg_path / "resources/Drug Screens/CTRP2_Res_nAUC.rds")[None]
        return bk_response.columns.to_list()
        

def spearman_r_pvals_vectorized(X, y, eps=1e-12):
    """Compute vectorized Spearman correlation p-values for all features.

    Args:
        X (array-like): Feature matrix with shape ``(n_samples, n_features)``.
        y (array-like): Response vector with shape ``(n_samples,)``.
        eps (float): Small constant used to avoid division by zero.

    Returns:
        np.ndarray: Two-sided p-values with shape ``(n_features,)``.
    """
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).reshape(-1)

    n = X.shape[0]
    if y.shape[0] != n:
        raise ValueError("X and y must have the same number of samples")

    # Rank-transform each feature and the response; rankdata handles ties.
    Xr = np.apply_along_axis(rankdata, 0, X)
    yr = rankdata(y)

    # Pearson correlation on ranks is Spearman correlation.
    Xc = Xr - Xr.mean(axis=0, keepdims=True)
    yc = yr - yr.mean()

    denom = (np.linalg.norm(Xc, axis=0) + eps) * (np.linalg.norm(yc) + eps)
    rho = (Xc.T @ yc) / denom
    rho = np.clip(rho, -1.0, 1.0)

    # Approximate the two-sided p-value using the t distribution.
    df = max(n - 2, 1)
    t = rho * np.sqrt(df / np.maximum(1.0 - rho**2, eps))
    pvals = 2.0 * t_dist.sf(np.abs(t), df=df)

    return pvals


def fdr_adjust(pvals, method="BH"):
    """Adjust p-values with BH or BY false-discovery-rate control.

    Args:
        pvals (array-like): Raw p-values.
        method (str): Adjustment method, either ``"BH"`` or ``"BY"``.

    Returns:
        np.ndarray: Adjusted q-values in the original input order.
    """
    p = np.asarray(pvals, dtype=np.float64)
    m = p.size
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, m + 1)

    if method.upper() == "BH":
        q = p * m / ranks
    elif method.upper() == "BY":
        c_m = np.sum(1.0 / np.arange(1, m + 1))
        q = p * m * c_m / ranks
    else:
        raise ValueError("method must be 'BH' or 'BY'")

    # Enforce monotonicity after sorting by raw p-value.
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    qvals = np.empty_like(q_sorted)
    qvals[order] = np.clip(q_sorted, 0.0, 1.0)
    return qvals
