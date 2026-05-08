"""
ensemble.py — Deep Ensemble architecture and training for Module 06.
深度集成模型架构与训练逻辑（模块 06 核心代码）。

Reference / 参考文献:
  Lakshminarayanan et al. (2017) "Simple and Scalable Predictive Uncertainty
  Estimation Using Deep Ensembles." NeurIPS 2017. arXiv:1612.01474.
"""

from __future__ import annotations

import logging
import random
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# EnsembleMember — single MLP with Gaussian output head
# 单个 ensemble 成员：3 层 MLP，输出高斯分布的均值和对数标准差
# ---------------------------------------------------------------------------

class EnsembleMember(nn.Module):
    """
    3-layer MLP that predicts a Gaussian distribution over a scalar binding score.
    输出一对 (mean, log_sigma)，训练时用 Gaussian NLL loss。

    hidden_dim: width of each hidden layer (default 256 per AGENT_TODO §3.5).
    dropout_p:  small dropout for regularisation; NOT used for MC Dropout inference.
    """

    def __init__(self, input_dim: int, hidden_dim: int = 256, dropout_p: float = 0.1):
        super().__init__()
        # 3-layer MLP backbone / 3 层 MLP 主干
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
        )
        # Separate output heads for mean and log_sigma
        # 分离的输出头：均值 + 对数标准差
        self.head_mean = nn.Linear(hidden_dim // 2, 1)
        self.head_log_sigma = nn.Linear(hidden_dim // 2, 1)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Returns (mean, sigma) — both shape (batch,).
        返回 (均值, 标准差)，各形状为 (batch,)。
        """
        h = self.net(x)
        mean = self.head_mean(h).squeeze(-1)
        # Clamp log_sigma to [-7, 7] for numerical stability
        # 将 log_sigma 截断到 [-7, 7] 防止数值爆炸（Lakshminarayanan 2017 §3.1）
        log_sigma = self.head_log_sigma(h).squeeze(-1).clamp(-7.0, 7.0)
        sigma = torch.exp(log_sigma)
        return mean, sigma


def gaussian_nll_loss(
    mean: torch.Tensor, sigma: torch.Tensor, target: torch.Tensor
) -> torch.Tensor:
    """
    Gaussian Negative Log-Likelihood loss.
    高斯负对数似然损失 = 0.5*(log(2π) + 2*log(σ) + ((y - μ)/σ)²)

    Using PyTorch's built-in for numerical stability; equivalent to the
    Lakshminarayanan 2017 formulation.
    """
    return nn.GaussianNLLLoss()(mean, target, sigma ** 2)


# ---------------------------------------------------------------------------
# Training utilities / 训练工具函数
# ---------------------------------------------------------------------------

def _set_seed(seed: int) -> None:
    """Set all random seeds for full reproducibility. 设置全部随机种子。"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def train_member(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
    seed: int,
    hidden_dim: int = 256,
    lr: float = 1e-4,
    max_epochs: int = 200,
    patience: int = 10,
    batch_size: int = 32,
) -> Tuple[EnsembleMember, Dict]:
    """
    Train a single EnsembleMember with early stopping on val NLL.
    用早停（val NLL）训练单个 ensemble 成员。

    Args:
        X_train, y_train: training features/labels (numpy arrays).
        X_val,   y_val:   validation features/labels.
        seed:  per-member random seed — controls init AND DataLoader shuffle.
        hidden_dim, lr, max_epochs, patience, batch_size: hyperparams.

    Returns:
        (trained model, training log dict)
        返回 (训练好的模型, 训练日志字典)
    """
    # Seed BEFORE model instantiation so weight init is reproducible per member
    # 在构建模型前设置种子，确保每个成员权重初始化不同
    _set_seed(seed)

    input_dim = X_train.shape[1]
    model = EnsembleMember(input_dim=input_dim, hidden_dim=hidden_dim)

    # Seed again before DataLoader so shuffle order also differs
    # 再次设置种子，确保 DataLoader shuffle 顺序也不同
    _set_seed(seed)

    X_t = torch.tensor(X_train, dtype=torch.float32)
    y_t = torch.tensor(y_train, dtype=torch.float32)
    X_v = torch.tensor(X_val, dtype=torch.float32)
    y_v = torch.tensor(y_val, dtype=torch.float32)

    train_ds = TensorDataset(X_t, y_t)
    # Use a seeded generator for DataLoader shuffle / DataLoader 用确定性 shuffle
    gen = torch.Generator()
    gen.manual_seed(seed)
    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True, generator=gen
    )

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best_val_loss = float("inf")
    best_state = None
    patience_counter = 0
    log: Dict = {
        "seed": seed,
        "train_nll": [],
        "val_nll": [],
        "best_epoch": 0,
        "final_val_nll": None,
    }

    for epoch in range(max_epochs):
        # --- Training pass / 训练 ---
        model.train()
        epoch_loss = 0.0
        for xb, yb in train_loader:
            optimizer.zero_grad()
            mean, sigma = model(xb)
            loss = gaussian_nll_loss(mean, sigma, yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(xb)
        train_nll = epoch_loss / len(X_train)

        # --- Validation pass / 验证 ---
        model.eval()
        with torch.no_grad():
            val_mean, val_sigma = model(X_v)
            val_nll = gaussian_nll_loss(val_mean, val_sigma, y_v).item()

        log["train_nll"].append(round(train_nll, 6))
        log["val_nll"].append(round(val_nll, 6))

        # Early stopping / 早停
        if val_nll < best_val_loss - 1e-6:
            best_val_loss = val_nll
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
            log["best_epoch"] = epoch
        else:
            patience_counter += 1
            if patience_counter >= patience:
                logger.info(
                    "Member seed=%d early stop at epoch %d (best val NLL=%.4f)",
                    seed, epoch, best_val_loss,
                )
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    log["final_val_nll"] = best_val_loss
    return model, log


# ---------------------------------------------------------------------------
# DeepEnsemble — container for N trained members
# 容器：管理 N 个训练好的成员
# ---------------------------------------------------------------------------

class DeepEnsemble:
    """
    Wraps N independently-trained EnsembleMember MLPs.
    每个成员用不同 seed 独立训练，保证多样性（diversity）。
    """

    def __init__(self, input_dim: int, hidden_dim: int = 256, n_members: int = 5):
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.n_members = n_members
        self.members: List[EnsembleMember] = []
        self.logs: List[Dict] = []
        self.meta: Dict = {}

    def train(
        self,
        X_train: np.ndarray,
        y_train: np.ndarray,
        X_val: np.ndarray,
        y_val: np.ndarray,
        **train_kwargs,
    ) -> None:
        """
        Train n_members independently (seeds 0 .. n_members-1).
        用 seed 0..N-1 分别独立训练 N 个成员。
        """
        t0 = time.time()
        self.members = []
        self.logs = []
        for seed in range(self.n_members):
            logger.info("Training member %d/%d (seed=%d)…", seed + 1, self.n_members, seed)
            model, log = train_member(
                X_train, y_train, X_val, y_val,
                seed=seed,
                hidden_dim=self.hidden_dim,
                **train_kwargs,
            )
            model.eval()
            self.members.append(model)
            self.logs.append(log)
        self.meta["train_time_s"] = round(time.time() - t0, 1)

    def predict(
        self, X: np.ndarray, return_epistemic: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Ensemble predictive mean and total uncertainty std.
        返回集成预测均值和总不确定性标准差。

        Total predictive variance = Var(means) + E[sigma²]
        总预测方差 = 各成员均值的方差（epistemic）+ 各成员 sigma² 的期望（aleatoric）
        This is the correct decomposition from Lakshminarayanan 2017 eq. (3).
        这是 Lakshminarayanan 2017 公式 (3) 的正确分解。

        Args:
            X:                Input array (N, input_dim).
            return_epistemic: If True, also return epistemic_std separately.
                              如果 True，额外返回纯 epistemic std。

        Returns:
            mean:       (N,) — ensemble mean.
            total_std:  (N,) — sqrt(epistemic_var + aleatoric_var).
            [epistemic_std]: (N,) — only if return_epistemic=True.
        """
        if not self.members:
            raise RuntimeError("Call train() before predict().")

        X_t = torch.tensor(X, dtype=torch.float32)
        all_means = []
        all_sigma2 = []
        for model in self.members:
            model.eval()
            with torch.no_grad():
                m, s = model(X_t)
                all_means.append(m.numpy())
                all_sigma2.append((s ** 2).numpy())

        stacked_means = np.stack(all_means, axis=0)   # (n_members, N)
        stacked_var = np.stack(all_sigma2, axis=0)    # (n_members, N)

        mean = stacked_means.mean(axis=0)
        epistemic_var = stacked_means.var(axis=0)     # Var of means
        aleatoric_var = stacked_var.mean(axis=0)      # E[sigma²]

        total_std = np.sqrt(np.clip(epistemic_var + aleatoric_var, 0.0, None))
        epistemic_std = np.sqrt(np.clip(epistemic_var, 0.0, None))

        if return_epistemic:
            return mean, total_std, epistemic_std  # type: ignore[return-value]
        return mean, total_std
