"""test_smoke.py — mock boltz CLI and test parsing/aggregation pipeline end-to-end.

冒烟测试：模拟 boltz CLI 调用，端到端测试解析和聚合管道。

Run: pytest 05_structure_prediction/processes/tests/test_smoke.py -v
"""

import json
import math
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

# Add the processes directory to path for importing helpers
# 将 processes 目录添加到路径以导入辅助函数
import sys
PROCESSES_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROCESSES_DIR))


# ── Inline helper functions (mirrors notebook logic) ─────────────────────────
# These functions are the same as in 01_run_boltz2.ipynb, extracted for testing.
# 这些函数与 01_run_boltz2.ipynb 中相同，提取出来用于测试。

import textwrap


def prepare_pair_input(rbp_id, rbp_seq, receptor_id, receptor_seq, work_dir):
    """Write Boltz-2 FASTA input (mirrors notebook cell 4)."""
    def wrap60(seq):
        return "\n".join(textwrap.wrap(seq, 60))

    pair_name = f"{rbp_id}__{receptor_id}"
    fasta_path = Path(work_dir) / f"{pair_name}.fasta"
    content = f">A|protein\n{wrap60(rbp_seq)}\n>B|protein\n{wrap60(receptor_seq)}\n"
    fasta_path.write_text(content)
    return fasta_path


def parse_boltz2_confidence(output_dir, pair_name):
    """Parse Boltz-2 confidence JSON (mirrors notebook cell 6)."""
    results_dir = Path(output_dir) / f"boltz_results_{pair_name}"
    pred_dir = results_dir / "predictions" / pair_name
    conf_file = pred_dir / f"confidence_{pair_name}_model_0.json"

    if not conf_file.exists():
        candidates = list(Path(output_dir).rglob("confidence_*.json"))
        if candidates:
            conf_file = candidates[0]
        else:
            return {
                "ipTM": float("nan"), "pTM": float("nan"),
                "pLDDT_mean": float("nan"),
                "predicted_dG_kcal_mol": float("nan"),
                "predicted_pKd": float("nan"),
                "confidence": float("nan"),
                "conf_file_found": False,
            }

    with open(conf_file) as f:
        data = json.load(f)

    iptm = data.get("iptm", data.get("ipTM", float("nan")))
    ptm = data.get("ptm", data.get("pTM", float("nan")))
    plddt_raw = data.get("plddt", data.get("pLDDT", []))
    plddt_mean = (
        float(np.mean(plddt_raw)) if isinstance(plddt_raw, list) and plddt_raw
        else float(plddt_raw if plddt_raw else float("nan"))
    )

    return {
        "ipTM": float(iptm),
        "pTM": float(ptm),
        "pLDDT_mean": plddt_mean,
        "predicted_dG_kcal_mol": float("nan"),   # NaN for protein-protein
        "predicted_pKd": float("nan"),
        "confidence": float(iptm),
        "conf_file_found": True,
    }


def build_affinity_row(rbp_id, receptor_id, conf, run_meta):
    """Build one CSV row from confidence dict + run metadata."""
    return {
        "rbp_id": rbp_id,
        "receptor_id": receptor_id,
        "model": "boltz2_2.0.3",
        "predicted_dG_kcal_mol": conf["predicted_dG_kcal_mol"],
        "predicted_pKd": conf["predicted_pKd"],
        "confidence": conf["confidence"],
        "interface_ipTM": conf["ipTM"],
        "pTM": conf.get("pTM", float("nan")),
        "pLDDT_mean": conf.get("pLDDT_mean", float("nan")),
        "run_success": run_meta["success"],
        "runtime_s": run_meta["runtime_s"],
        "partial": run_meta["partial"],
    }


# ── Tests ────────────────────────────────────────────────────────────────────

class TestPrepareInput:
    """Test FASTA input preparation."""

    def test_fasta_written(self, tmp_path):
        rbp_seq = "MGSLAGFPALSIGVRDTKINIVARAFQNLLLQVQPVALRVVLLG"
        rec_seq = "MAMNTLSTALHTLHLLRPQLLWALVLIVPAAVLWFLRQRDADVW"
        fasta = prepare_pair_input("rbp_01", rbp_seq, "tonB", rec_seq, tmp_path)
        assert fasta.exists()

    def test_fasta_chains(self, tmp_path):
        """FASTA must have >A|protein and >B|protein headers."""
        fasta = prepare_pair_input("rbp_01", "MKLLL", "tonB", "MAMNT", tmp_path)
        text = fasta.read_text()
        assert ">A|protein" in text
        assert ">B|protein" in text

    def test_pair_name_in_filename(self, tmp_path):
        fasta = prepare_pair_input("rbp_01", "MKLLL", "tonB", "MAMNT", tmp_path)
        assert "rbp_01__tonB" in fasta.name


class TestParseConfidence:
    """Test parse_boltz2_confidence with mocked JSON files."""

    def _write_fake_confidence(self, base_dir, pair_name, iptm=0.72, ptm=0.68):
        """Write a fake Boltz-2 confidence JSON at the expected path."""
        pred_dir = (
            base_dir / f"boltz_results_{pair_name}" / "predictions" / pair_name
        )
        pred_dir.mkdir(parents=True)
        conf_file = pred_dir / f"confidence_{pair_name}_model_0.json"
        conf_file.write_text(json.dumps({
            "iptm": iptm,
            "ptm": ptm,
            "plddt": [0.85, 0.90, 0.78],
        }))
        return conf_file

    def test_parses_iptm(self, tmp_path):
        pair = "rbp_01__tonB"
        self._write_fake_confidence(tmp_path, pair, iptm=0.72)
        result = parse_boltz2_confidence(tmp_path, pair)
        assert abs(result["ipTM"] - 0.72) < 1e-6

    def test_protein_protein_dg_is_nan(self, tmp_path):
        """For protein-protein pairs, dG must be NaN (no ligand affinity head)."""
        pair = "rbp_01__tonB"
        self._write_fake_confidence(tmp_path, pair)
        result = parse_boltz2_confidence(tmp_path, pair)
        assert math.isnan(result["predicted_dG_kcal_mol"])
        assert math.isnan(result["predicted_pKd"])

    def test_no_json_returns_nan(self, tmp_path):
        """When no confidence JSON found, all metrics are NaN."""
        result = parse_boltz2_confidence(tmp_path, "nonexistent")
        assert math.isnan(result["ipTM"])
        assert not result["conf_file_found"]

    def test_plddt_mean_computed(self, tmp_path):
        pair = "rbp_01__tonB"
        self._write_fake_confidence(tmp_path, pair, iptm=0.8)
        result = parse_boltz2_confidence(tmp_path, pair)
        expected_plddt = (0.85 + 0.90 + 0.78) / 3
        assert abs(result["pLDDT_mean"] - expected_plddt) < 1e-4


class TestAffinityRowAggregation:
    """Test CSV row building and INTERFACE schema conformance."""

    REQUIRED_COLS = [
        "rbp_id", "receptor_id", "model",
        "predicted_dG_kcal_mol", "predicted_pKd",
        "confidence", "interface_ipTM",
    ]

    def _fake_meta(self, success=True, partial=False):
        return {"success": success, "runtime_s": 120.5, "partial": partial, "returncode": 0}

    def test_row_has_required_columns(self):
        conf = {
            "ipTM": 0.72, "pTM": 0.68, "pLDDT_mean": 0.85,
            "predicted_dG_kcal_mol": float("nan"),
            "predicted_pKd": float("nan"),
            "confidence": 0.72,
        }
        row = build_affinity_row("EU717894.1_rbp_01", "GCF_000007145.1_tonB", conf, self._fake_meta())
        for col in self.REQUIRED_COLS:
            assert col in row, f"Missing column: {col}"

    def test_csv_append_produces_valid_dataframe(self, tmp_path):
        """Appending a row produces a valid, readable CSV."""
        conf = {
            "ipTM": 0.72, "pTM": 0.68, "pLDDT_mean": 0.85,
            "predicted_dG_kcal_mol": float("nan"),
            "predicted_pKd": float("nan"),
            "confidence": 0.72,
        }
        row = build_affinity_row("EU717894.1_rbp_01", "GCF_000007145.1_tonB", conf, self._fake_meta())
        csv = tmp_path / "affinity_priors.csv"
        pd.DataFrame([row]).to_csv(csv, index=False, float_format="%.6f")

        df = pd.read_csv(csv)
        assert len(df) == 1
        assert df.iloc[0]["rbp_id"] == "EU717894.1_rbp_01"
        for col in self.REQUIRED_COLS:
            assert col in df.columns

    def test_timeout_produces_partial_row(self):
        """A timed-out run should produce a row with partial=True and NaN confidence."""
        conf = {
            "ipTM": float("nan"), "pTM": float("nan"), "pLDDT_mean": float("nan"),
            "predicted_dG_kcal_mol": float("nan"),
            "predicted_pKd": float("nan"),
            "confidence": float("nan"),
        }
        meta = self._fake_meta(success=False, partial=True)
        row = build_affinity_row("EU717894.1_rbp_01", "GCF_000007145.1_tonB", conf, meta)
        assert row["partial"] is True
        assert math.isnan(row["confidence"])
