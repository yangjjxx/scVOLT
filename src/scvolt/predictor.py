"""ElasticNet-based predictor for scVOLT latent factors."""

import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import ElasticNetCV


class ElasticNetPredictor:
    """Fit and apply an ElasticNet regression model for sensitivity scores."""

    def __init__(self, max_iter: int = 5000, cv: int = 5, n_jobs: int = 10):
        """Store ElasticNet cross-validation settings.

        Args:
            max_iter (int): Maximum solver iterations for ElasticNetCV.
            cv (int): Number of cross-validation folds.
            n_jobs (int): Number of parallel jobs used by ElasticNetCV.
        """
        self.max_iter = max_iter
        self.cv = cv
        self.n_jobs = n_jobs
        self.reg = None

    def train(self, X, y):
        """Train the ElasticNet model on latent features and response values.

        Args:
            X (array-like): Training feature matrix.
            y (array-like): Training response vector.
        """
        X = np.asarray(X)
        y = np.asarray(y).reshape(-1)
        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y must have the same number of samples")
        if y.shape[0] < 2:
            raise ValueError("ElasticNet training requires at least two samples")

        # Keep CV valid for small datasets while respecting the configured upper bound.
        cv = min(self.cv, y.shape[0])

        # Build the regression pipeline before fitting so CV settings stay local to this call.
        self.reg = make_pipeline(
            # StandardScaler can be enabled if future inputs are not already normalized.
            ElasticNetCV(
                l1_ratio=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                alphas=[1e-4, 1e-3, 1e-2, 1e-1, 1, 10],
                fit_intercept=True,
                max_iter=self.max_iter,
                cv=cv,
                n_jobs=self.n_jobs,
            )
        )
        self.reg.fit(X, y)

    def predict(self, X_test) -> np.ndarray:
        """Predict clipped sensitivity scores for new latent features.

        Args:
            X_test (array-like): Feature matrix to score.

        Returns:
            np.ndarray: Predicted scores clipped to the ``[0, 1]`` range.
        """
        # Clip predictions to the response scale used by the drug screens.
        y_pred = self.reg.predict(X_test)
        return np.clip(y_pred, 0, 1)

    def best_params(self):
        """Return the selected ElasticNet hyperparameters after training.

        Returns:
            dict: The selected alpha and l1_ratio values.
        """
        model = self.reg.named_steps["elasticnetcv"]
        return {
            "alpha": model.alpha_,
            "l1_ratio": model.l1_ratio_,
        }
