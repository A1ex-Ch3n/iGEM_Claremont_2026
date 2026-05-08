"""test_schema.py — verify that module 05 output files match INTERFACE.md §Module 05.

模式测试：验证模块05的输出文件符合 INTERFACE.md §Module 05 的列规范。

Run: pytest 05_structure_prediction/processes/tests/test_schema.py -v
"""

import json
import math
from pathlib import Path

import pandas as pd
import pytest

# Anchor to module root: tests/ → processes/ → 05_structure_prediction/
# 路径锚点：tests/ → processes/ → 05_structure_prediction/
MODULE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT   = MODULE_ROOT.parent
OUTPUTS_DIR = MODULE_ROOT / "outputs"
BOLTZ_DIR = OUTPUTS_DIR / "boltz2"


# ── INTERFACE §Module 05: required columns for affinity_priors.csv ──────────
REQUIRED_AFFINITY_COLS = [
    "rbp_id",
    "receptor_id",
    "model",
    "predicted_dG_kcal_mol",
    "predicted_pKd",
    "confidence",
    "interface_ipTM",
]

# Required keys in a Boltz-2 affinity JSON (when present)
# Boltz-2 亲和力 JSON 的必需键
REQUIRED_AFFINITY_JSON_KEYS = [
    "predicted_dG_kcal_mol",
    "predicted_pKd",
    "confidence",
]


class TestAffinityPriorsSchema:
    """Verify affinity_priors.csv matches INTERFACE.md column contract."""

    def test_file_exists(self):
        """affinity_priors.csv must exist after any run."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        assert csv.exists(), (
            f"affinity_priors.csv missing at {csv}. "
            "Run 01_run_boltz2.ipynb to generate."
        )

    def test_required_columns(self):
        """All INTERFACE-required columns must be present."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        missing = [c for c in REQUIRED_AFFINITY_COLS if c not in df.columns]
        assert not missing, f"Missing columns: {missing}"

    def test_no_empty_ids(self):
        """rbp_id and receptor_id must not be empty strings."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        assert df["rbp_id"].notna().all(), "NaN in rbp_id"
        assert df["receptor_id"].notna().all(), "NaN in receptor_id"
        assert (df["rbp_id"] != "").all(), "Empty string in rbp_id"

    def test_confidence_in_range(self):
        """confidence values must be NaN or in [0, 1]."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        valid = df["confidence"].apply(
            lambda x: math.isnan(x) or 0.0 <= x <= 1.0
        )
        assert valid.all(), (
            f"confidence out of [0,1]: {df.loc[~valid, 'confidence'].tolist()}"
        )

    def test_interface_iptm_in_range(self):
        """interface_ipTM must be NaN or in [0, 1]."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        valid = df["interface_ipTM"].apply(
            lambda x: math.isnan(x) or 0.0 <= x <= 1.0
        )
        assert valid.all(), "interface_ipTM out of [0, 1]"

    def test_model_string_not_empty(self):
        """model column must be a non-empty string."""
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        assert (df["model"].str.len() > 0).all(), "Empty model string"


class TestBoltz2AffinityJson:
    """Verify the affinity.json written by the notebook matches the contract.

    The notebook writes a custom affinity.json; this test validates it.
    Notebook 写入自定义 affinity.json，本测试验证其格式。
    """

    def test_affinity_json_keys(self):
        """Each Boltz-2 output dir should contain an affinity.json with required keys."""
        aff_jsons = list(BOLTZ_DIR.rglob("affinity.json"))
        if not aff_jsons:
            pytest.skip("No affinity.json found — Boltz-2 not yet run")
        for jf in aff_jsons:
            data = json.loads(jf.read_text())
            missing = [k for k in REQUIRED_AFFINITY_JSON_KEYS if k not in data]
            assert not missing, f"{jf}: missing keys {missing}"


class TestManifestSchema:
    """Verify outputs/MANIFEST.csv matches INTERFACE.md universal convention."""

    REQUIRED_MANIFEST_COLS = [
        "filename",
        "sha256",
        "bytes",
        "n_records",
        "created_utc",
        "source_acc",
        "source_module",
        "notes",
    ]

    def test_manifest_exists(self):
        manifest = OUTPUTS_DIR / "MANIFEST.csv"
        assert manifest.exists(), (
            f"MANIFEST.csv missing at {manifest}. "
            "Run 01_run_boltz2.ipynb to generate."
        )

    def test_manifest_columns(self):
        manifest = OUTPUTS_DIR / "MANIFEST.csv"
        if not manifest.exists():
            pytest.skip("MANIFEST.csv not yet generated")
        df = pd.read_csv(manifest)
        missing = [c for c in self.REQUIRED_MANIFEST_COLS if c not in df.columns]
        assert not missing, f"MANIFEST.csv missing columns: {missing}"
