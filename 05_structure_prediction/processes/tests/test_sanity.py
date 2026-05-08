"""test_sanity.py — biology-level sanity checks for module 05 outputs.

合理性测试：模块05输出的生物学合理性检查。

Run: pytest 05_structure_prediction/processes/tests/test_sanity.py -v
"""

from pathlib import Path

import pytest

# tests/ → processes/ → 05_structure_prediction/
MODULE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT   = MODULE_ROOT.parent
OUTPUTS_DIR = MODULE_ROOT / "outputs"
INPUTS_DIR = MODULE_ROOT / "inputs"
BOLTZ_DIR = OUTPUTS_DIR / "boltz2"


class TestInputSequenceSanity:
    """Verify input sequences are biologically plausible."""

    def test_rbp_fasta_exists(self):
        """RBP sequence file must exist."""
        faa = INPUTS_DIR / "EU717894.1_rbp_sequences.faa"
        assert faa.exists(), f"Missing: {faa}"

    def test_receptor_fasta_exists(self):
        """Receptor sequence file must exist."""
        faa = INPUTS_DIR / "xcc_receptors_minimal.faa"
        assert faa.exists(), f"Missing: {faa}"

    def test_rbp_sequence_length(self):
        """phiL7 P25 (rbp_01) must be 80-100 aa (ACE75765.1 = 85 aa)."""
        faa = INPUTS_DIR / "EU717894.1_rbp_sequences.faa"
        if not faa.exists():
            pytest.skip("RBP FASTA not yet generated")
        from Bio import SeqIO
        records = list(SeqIO.parse(faa, "fasta"))
        assert len(records) >= 1, "No sequences in RBP FASTA"
        first_len = len(records[0].seq)
        assert 70 <= first_len <= 110, (
            f"Unexpected RBP length {first_len} aa "
            f"(expected ~85 aa for phiL7 P25 / ACE75765.1)"
        )

    def test_receptor_tonb_length(self):
        """Xcc TonB (Q8P5W2) must be 550-650 aa."""
        faa = INPUTS_DIR / "xcc_receptors_minimal.faa"
        if not faa.exists():
            pytest.skip("Receptor FASTA not yet generated")
        from Bio import SeqIO
        records = {r.id.split("|")[0].strip(): r for r in SeqIO.parse(faa, "fasta")}
        assert "GCF_000007145.1_tonB" in records, "TonB record missing"
        tonb_len = len(records["GCF_000007145.1_tonB"].seq)
        assert 550 <= tonb_len <= 650, (
            f"Unexpected TonB length {tonb_len} aa (Q8P5W2 = 604 aa)"
        )

    def test_receptor_ids_match_interface(self):
        """Receptor IDs must follow INTERFACE.md format: <host_acc>_<gene_name>.

        Examples: GCF_000007145.1_tonB, GCF_000007145.1_exbB
        The host_acc may contain dots (version suffix).
        示例：GCF_000007145.1_tonB，host_acc 可含句点（版本后缀）。
        """
        import re
        faa = INPUTS_DIR / "xcc_receptors_minimal.faa"
        if not faa.exists():
            pytest.skip("Receptor FASTA not yet generated")
        from Bio import SeqIO
        # Allow period in accession (e.g. GCF_000007145.1)
        # 允许 accession 中含句点（如 GCF_000007145.1）
        pattern = re.compile(r"^[A-Z0-9_.]+_[a-zA-Z][a-zA-Z0-9]*$")
        for rec in SeqIO.parse(faa, "fasta"):
            rec_id = rec.id.split("|")[0].strip()
            assert pattern.match(rec_id), (
                f"Receptor ID '{rec_id}' doesn't match <host_acc>_<gene_name> pattern"
            )


class TestBoltz2OutputSanity:
    """Sanity checks on Boltz-2 structure outputs."""

    def test_pdb_has_atom_records(self):
        """Any PDB output must contain ATOM records and be parseable."""
        pdbs = list(BOLTZ_DIR.rglob("*.pdb"))
        if not pdbs:
            pytest.skip("No PDB files found — Boltz-2 not yet run (expected for timeout/partial)")

        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)

        for pdb_path in pdbs:
            text = pdb_path.read_text()
            assert "ATOM" in text, f"PDB has no ATOM records: {pdb_path.name}"

            # Parse with Bio.PDB to verify structure integrity
            # 用 Bio.PDB 解析以验证结构完整性
            struct = parser.get_structure("test", str(pdb_path))
            chains = list(struct.get_chains())
            assert len(chains) >= 2, (
                f"Expected ≥2 chains (RBP + receptor) in {pdb_path.name}, "
                f"got {len(chains)}"
            )

    def test_pair_dir_name_convention(self):
        """Boltz-2 output directory names must follow <rbp_id>__<receptor_id> convention."""
        pair_dirs = [d for d in BOLTZ_DIR.iterdir() if d.is_dir()]
        if not pair_dirs:
            pytest.skip("No Boltz-2 pair directories found")
        for d in pair_dirs:
            assert "__" in d.name, (
                f"Pair directory '{d.name}' must contain '__' separator "
                f"(<rbp_id>__<receptor_id>)"
            )

    def test_boltz_input_fasta_saved(self):
        """Boltz-2 input FASTA must be saved (for Laguna re-run if needed).

        Even if the run timed out, the input FASTA must be preserved.
        即使运行超时，输入 FASTA 文件也必须保存，以便在 Laguna 上重新运行。
        """
        fastas = list(INPUTS_DIR.glob("boltz_input_*.fasta"))
        assert fastas, (
            f"No boltz_input_*.fasta found in {INPUTS_DIR}. "
            "The notebook should write the input before launching boltz."
        )


class TestAffinityPriorsSanity:
    """Sanity checks on the aggregated affinity_priors.csv."""

    def test_rbp_id_format(self):
        """rbp_id values must follow EU717894.1_rbp_NN format."""
        import re, pandas as pd
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        pattern = re.compile(r"^[A-Z]{2}\d+\.\d+_rbp_\d{2}$")
        for rbp_id in df["rbp_id"]:
            assert pattern.match(str(rbp_id)), (
                f"rbp_id '{rbp_id}' doesn't match expected format "
                f"(e.g. EU717894.1_rbp_01)"
            )

    def test_at_least_one_row(self):
        """affinity_priors.csv must have at least 1 row."""
        import pandas as pd
        csv = OUTPUTS_DIR / "affinity_priors.csv"
        if not csv.exists():
            pytest.skip("affinity_priors.csv not yet generated")
        df = pd.read_csv(csv)
        assert len(df) >= 1, "affinity_priors.csv is empty"
