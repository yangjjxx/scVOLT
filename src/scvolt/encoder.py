"""Matrix factorization encoder with entropic optimal transport losses."""

from typing import Callable, List
import numpy as np
import torch
import torch.nn.functional as F
from scvolt import utils
from torch import optim
from tqdm import tqdm


class MatrixFactorizationEncoder:
    """Perform integrative non-negative matrix factorization with OT regularization.

    Args:
        latent_dim (int): Latent dimension for the shared and resolution-specific factors.
        h_regularization (dict | float): Entropy regularization for each H factor.
        w_regularization (float): Entropy regularization for the shared W factor.
        eps (float): Entropic regularization for optimal transport.
        cost (str): SciPy distance metric used for the empirical ground cost.
    """

    def __init__(
            self,
            latent_dim: int,
            h_regularization: dict,
            w_regularization: float,
            eps: float,
            cost: str
    ):
        """Store encoder hyperparameters and initialize training state.

        Args:
            latent_dim (int): Latent dimension for factorization.
            h_regularization (dict | float): Entropy regularization for H factors.
            w_regularization (float): Entropy regularization for W.
            eps (float): Entropic regularization for OT.
            cost (str): Ground-cost metric passed to SciPy ``cdist``.
        """

        # Save hyperparameters for later optimization steps.
        self.latent_dim = latent_dim
        self.w_regularization = w_regularization
        self.h_regularization = h_regularization
        self.eps = eps
        self.cost = cost

        # Track loss histories across alternating optimization steps.
        self.losses_w, self.losses_h, self.losses = [], [], []

        # Store per-resolution input matrices, factors, dual variables, and kernels.
        self.A, self.H, self.G, self.K = {}, {}, {}, {}

    def init_parameters(
            self,
            rna_data: dict,
            dtype: torch.dtype,
            device: torch.device,
            normalize_rows: bool = False,
            enable_gene_weight: bool = False,
            gene_weight: dict = None
    ) -> None:
        """Initialize factors, dual variables, and transport kernels from input data.

        Args:
            rna_data (dict): Mapping from resolution name to expression matrix.
            dtype (torch.dtype): Tensor dtype used by PyTorch.
            device (torch.device): Tensor device used by PyTorch.
            normalize_rows (bool): Reserved flag; current implementation normalizes columns.
            enable_gene_weight (bool): Whether to use per-gene weights.
            gene_weight (dict): Optional per-resolution gene weights.
        """

        # Record dimensions and regularization constants shared by all resolutions.
        self.resolutions = rna_data.keys()
        self.n_resolutions = len(rna_data)
        self.n_obs = {}
        self.n_vars = len(rna_data["bk"].columns)
        self.mu = {}
        self.rho = self.w_regularization / (np.log(self.n_vars) * self.latent_dim)

        # Initialize the per-resolution gene-weight dictionary.
        self.gene_weight = {}

        if not isinstance(self.h_regularization, dict):
            self.h_regularization = {resolution: self.h_regularization for resolution in self.resolutions}

        # Initialize per-resolution matrices and factors.
        for resolution in self.resolutions:
            # Convert the expression matrix into a feature-by-sample tensor.
            self.A[resolution] = utils.reference_dataset(rna_data[resolution].to_numpy(), dtype, device)
            self.n_obs[resolution] = self.A[resolution].shape[1]

            # Compute the resolution-specific entropy coefficient.
            self.mu[resolution] = self.h_regularization[resolution] / (np.log(self.n_obs[resolution]) * self.latent_dim)

            # Use supplied gene weights or default to equal weights.
            if enable_gene_weight:
                weight_matrix = gene_weight[resolution].to_numpy()
                weight_matrix = torch.Tensor(weight_matrix).reshape(-1, 1)
                weight_matrix = weight_matrix.to(dtype=dtype, device=device)
                self.gene_weight[resolution] = weight_matrix
            else:
                self.gene_weight[resolution] = torch.ones(
                    self.n_vars, 1, dtype=dtype, device=device
                )

            # Add a small value for numerical stability, then normalize each sample column.
            self.A[resolution] += 1e-6
            self.A[resolution] /= self.A[resolution].sum(0)

            # Initialize and normalize the resolution-specific factor H.
            self.H[resolution] = torch.ones(
                self.latent_dim, self.n_obs[resolution],  device=device, dtype=dtype
            )
            self.H[resolution] = utils.normalize_tensor(self.H[resolution])

            # Initialize the OT dual variable G for this resolution.
            self.G[resolution] = torch.zeros_like(self.A[resolution], requires_grad=True)
        features = None
        for resolution in self.resolutions:
            # Build the ground-cost kernel from normalized expression features.
            features = 1e-6+self.A[resolution].cpu().numpy()
            self.K[resolution] = utils.compute_ground_cost(
                features, self.cost, self.eps, dtype, device
            )

        # Initialize and normalize the shared factor W.
        self.W = torch.rand(self.n_vars, self.latent_dim, device=device, dtype=dtype)
        self.W = utils.normalize_tensor(self.W)

        # Release the temporary NumPy feature array.
        del features

    def train(
            self,
            rna_data: dict,
            max_iter_inner: int,
            max_iter: int,
            device: torch.device,
            dtype: torch.dtype,
            lr: float,
            optim_name: str,
            tol_inner: float,
            tol_outer: float,
            normalize_rows: bool,
            enable_gene_weight: bool,
            gene_weight: dict
    ) -> dict:
        """Train the encoder with alternating updates of H and W.

        Args:
            rna_data (dict): Input expression matrices by resolution.
            max_iter_inner (int): Maximum inner optimizer iterations for each factor update.
            max_iter (int): Maximum alternating update iterations.
            device (torch.device): Tensor device used by PyTorch.
            dtype (torch.dtype): Tensor dtype used by PyTorch.
            lr (float): Learning rate passed to the optimizer.
            optim_name (str): Optimizer name, one of ``"lbfgs"``, ``"sgd"``, or ``"adam"``.
            tol_inner (float): Inner-loop early-stopping tolerance.
            tol_outer (float): Outer-loop early-stopping tolerance.
            normalize_rows (bool): Reserved normalization flag passed to initialization.
            enable_gene_weight (bool): Whether to use per-gene weights.
            gene_weight (dict): Optional per-resolution gene weights.

        Returns:
            dict: Latent factors keyed by ``"H_<resolution>"`` plus shared ``"W"``.
        """

        # Initialize all factors, transport kernels, and dual variables.
        self.init_parameters(
            rna_data,
            dtype=dtype,
            device=device,
            normalize_rows=normalize_rows,
            enable_gene_weight=enable_gene_weight,
            gene_weight=gene_weight
        )

        self.lr = lr
        self.optim_name = optim_name

        # Reset loss histories for this training run.
        self.losses_h, self.losses_w, self.losses = [], [], []

        # Track one progress tick for each H and W update.
        pbar = tqdm(total=2 * max_iter, position=0, leave=True)

        # Alternate between updating the resolution-specific H factors and shared W.
        try:
            for _ in range(max_iter):

                # Optimize the dual variables for the H update.
                self.optimize(
                    loss_fn=self.loss_fn_h,
                    device=device,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    history=self.losses_h,
                    pbar=pbar,
                )

                # Recompute each resolution-specific H from the current dual variables.
                for resolution in self.resolutions:
                    wtg = self.W.T @ (self.G[resolution].detach() * self.gene_weight[resolution])
                    wtg /= self.n_obs[resolution]
                    # softmin incorporates the required negative sign.
                    wtg = wtg - wtg.mean(dim=0)
                    self.H[resolution] = F.softmin(wtg / self.mu[resolution], dim=0)

                pbar.update(1)

                # Save the outer-loop loss after updating H.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Optimize the dual variables for the W update.
                self.optimize(
                    loss_fn=self.loss_fn_w,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    history=self.losses_w,
                    pbar=pbar,
                    device=device,
                )

                # Aggregate all resolution-specific contributions before updating W.
                ght = 0
                for resolution in self.resolutions:
                    ght += (self.G[resolution] * self.gene_weight[resolution]) @ self.H[resolution].T/ self.n_obs[resolution]
                ght = ght - ght.mean(dim=0)
                self.W = F.softmin(ght.detach() / (self.n_resolutions * self.rho), dim=0)

                del ght

                pbar.update(1)

                # Save the outer-loop loss after updating W.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Stop when the monitored outer loss has converged.
                if utils.early_stop(self.losses, tol_outer, nonincreasing=True):
                    break

        except KeyboardInterrupt:
            print("Training interrupted.")

        # Return the learned factors as NumPy arrays for downstream predictors.
        ot_dict = {}
        for resolution in self.resolutions:
            ot_dict["H_" + resolution] = self.H[resolution].T.cpu().numpy()
            ot_dict["W"] = self.W.T.cpu().numpy()

        return ot_dict

    def build_optimizer(
            self, params, lr: float, optim_name: str
    ) -> torch.optim.Optimizer:
        """Create the requested PyTorch optimizer for the dual variables.

        Args:
            params (Iterable[torch.Tensor]): Parameters to optimize.
            lr (float): Optimizer learning rate.
            optim_name (str): Optimizer name, among ``"lbfgs"``, ``"sgd"``, and ``"adam"``.

        Returns:
            torch.optim.Optimizer: The optimizer.
        """
        # LBFGS is configured for one internal step because the caller controls the loop.
        if optim_name == "lbfgs":
            return optim.LBFGS(
                params,
                lr=lr,
                history_size=5,
                max_iter=1,
                line_search_fn="strong_wolfe",
            )
        elif optim_name == "sgd":
            return optim.SGD(params, lr=lr)
        elif optim_name == "adam":
            return optim.Adam(params, lr=lr)
        raise ValueError("optim_name must be one of 'lbfgs', 'sgd', or 'adam'")

    def optimize(
            self,
            loss_fn: Callable,
            max_iter: int,
            history: List,
            tol: float,
            pbar,
            device: torch.device,
    ) -> None:
        """Optimize a given function.

        Args:
            loss_fn (Callable): The function to optimize.
            max_iter (int): The maximum number of iterations.
            history (List): A list to append the losses to.
            tol (float): The tolerance before early stopping.
            pbar (A tqdm progress bar): The progress bar.
            device (str): The device to work on.
        """

        # Build an optimizer over all resolution-specific dual variables.
        optimizer = self.build_optimizer(
            [self.G[resolution] for resolution in self.G], lr=self.lr, optim_name=self.optim_name
        )

        # Display the most recent outer loss while the inner optimizer runs.
        if len(self.losses) > 0:
            total_loss = self.losses[-1].cpu().numpy()
        else:
            total_loss = "?"

        # Run the inner optimization loop until convergence or max_iter.
        for i in range(max_iter):
            def closure():
                """Evaluate the loss and gradients for one optimizer step."""
                optimizer.zero_grad()
                loss = loss_fn()
                loss.backward()
                return loss


            # Perform an optimizer step using the closure required by LBFGS.
            optimizer.step(closure)
            # Every ten steps, record the inner loss and refresh progress metadata.
            if i % 10 == 0:

                history.append(loss_fn().cpu().detach())

                pbar.set_postfix(
                    {
                        "loss": total_loss,
                        "loss_inner": history[-1].cpu().numpy(),
                        "inner_steps": i,
                    }
                )

                # Stop the inner loop once its loss has converged.
                if utils.early_stop(history, tol):
                    break

    @torch.no_grad()
    def total_dual_loss(self) -> torch.Tensor:
        """Compute the total dual loss for reporting and early stopping.

        Returns:
            torch.Tensor: The total dual loss.
        """

        # Accumulate each resolution's OT, constraint, and entropy terms.
        loss = 0

        for resolution in self.resolutions:
            # Subtract the normalized OT dual term.
            loss -= (
                    utils.ot_dual_loss(
                        self.A[resolution],
                        self.G[resolution],
                        self.K[resolution],
                        self.eps,
                        self.gene_weight[resolution],
                    )
                    / self.n_obs[resolution]
            )

            # Add the Lagrange multiplier term that links W and H.
            lagrange = self.W @ self.H[resolution]
            lagrange *= self.gene_weight[resolution] * self.G[resolution]
            lagrange = lagrange.sum()
            loss += lagrange / self.n_obs[resolution]

            # Subtract the entropy penalty for the current H factor.
            loss -= self.mu[resolution] * utils.entropy(self.H[resolution], min_one=True)

        # Subtract the entropy penalty for W.
        loss -= (
                self.n_resolutions * self.rho
                * utils.entropy(self.W, min_one=True)
        )
        return loss

    def loss_fn_h(self) -> torch.Tensor:
        """Compute the inner loss used while updating H.

        Returns:
            torch.Tensor: The H-update loss.
        """
        loss_h = 0
        for resolution in self.resolutions:
            # Add the normalized OT dual loss.
            odl = (
                    utils.ot_dual_loss(
                        self.A[resolution],
                        self.G[resolution],
                        self.K[resolution],
                        self.eps,
                        self.gene_weight[resolution],
                    )
                    / self.n_obs[resolution]
            )
            loss_h += odl

            # Subtract the entropy dual term induced by the current W factor.
            wtg = self.W.T @ (self.gene_weight[resolution] * self.G[resolution])
            wtg /= self.n_obs[resolution]
            edl = self.mu[resolution] * utils.entropy_dual_loss(-wtg / self.mu[resolution])
            loss_h -= edl
            del wtg

        return loss_h

    def loss_fn_w(self) -> torch.Tensor:
        """Compute the inner loss used while updating W.

        Returns:
            torch.Tensor: The W-update loss.
        """
        loss_w, ght = 0, 0

        for resolution in self.resolutions:

            # Add the normalized OT dual loss.
            odl = (
                    utils.ot_dual_loss(
                        self.A[resolution],
                        self.G[resolution],
                        self.K[resolution],
                        self.eps,
                        self.gene_weight[resolution],
                    )
                    / self.n_obs[resolution]
            )
            loss_w += odl

            # Accumulate the H-weighted dual variables for the entropy dual term.
            ght += (self.G[resolution] * self.gene_weight[resolution]) @ self.H[resolution].T / self.n_obs[resolution]

        # Subtract the shared W entropy dual term.
        edl = self.n_resolutions * self.rho * utils.entropy_dual_loss(-ght / (self.n_resolutions * self.rho))

        del ght

        loss_w -= edl

        return loss_w

