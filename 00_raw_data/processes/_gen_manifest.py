"""
Standalone script that implements the same logic as 01_verify_dataset.ipynb.
Called during notebook development to verify correctness before packaging.
此脚本实现与 01_verify_dataset.ipynb 相同的逻辑，用于 notebook 开发期间的正确性验证。
"""
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Path anchoring per INTERFACE.md §Universal conventions
# 路径锚定遵循 INTERFACE.md §Universal conventions
MODULE_ROOT = Path(__file__).resolve().parents[1]
PHAGE_ROOT = MODULE_ROOT / "phage"
BACTERIA_ROOT = MODULE_ROOT / "bacteria"


def sha256_file(fpath: Path, chunk_size: int = 65536) -> str:
    """Chunked SHA-256 computation. / 分块读取计算 SHA-256。"""
    h = hashlib.sha256()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def count_fasta_records(fpath: Path) -> float:
    """Count FASTA records; NaN if not a FASTA file. / 计算 FASTA 条数；非 FASTA 文件返回 NaN。"""
    if fpath.suffix.lower() not in (".fna", ".faa", ".fa", ".fasta"):
        return float("nan")
    try:
        return float(sum(1 for _ in SeqIO.parse(str(fpath), "fasta")))
    except Exception:
        return float("nan")


def main():
    print("=== Module 00 dataset verification / 模块 00 数据集验证 ===")

    # ------------------------------------------------------------------ #
    # Step 1: Walk file tree, build inventory                              #
    # 步骤 1：遍历文件树，建立清单                                          #
    # ------------------------------------------------------------------ #
    records = []
    for taxon_dir in [PHAGE_ROOT, BACTERIA_ROOT]:
        for acc_dir in sorted(taxon_dir.iterdir()):
            if not acc_dir.is_dir():
                continue
            source_acc = acc_dir.name
            for fpath in sorted(acc_dir.iterdir()):
                if fpath.is_file():
                    records.append(
                        {
                            "filename": str(fpath.relative_to(MODULE_ROOT)),
                            "bytes": fpath.stat().st_size,
                            "source_acc": source_acc,
                            "_path": fpath,
                        }
                    )

    inv_df = pd.DataFrame(records)
    print(f"Total files found: {len(inv_df)}")

    # ------------------------------------------------------------------ #
    # Step 2: Count FASTA records                                          #
    # 步骤 2：计算每个 FASTA 文件中的序列条数                               #
    # ------------------------------------------------------------------ #
    print("Counting FASTA records... / 计算序列条数...")
    inv_df["n_records"] = inv_df["_path"].apply(count_fasta_records)

    # ------------------------------------------------------------------ #
    # Step 3: Compute SHA-256                                              #
    # 步骤 3：计算 SHA-256 校验和                                           #
    # ------------------------------------------------------------------ #
    print("Computing SHA-256... / 计算 SHA-256（可能需要几分钟）...")
    inv_df["sha256"] = inv_df["_path"].apply(sha256_file)

    # ------------------------------------------------------------------ #
    # Step 4: Build MANIFEST.csv                                          #
    # 步骤 4：构建 MANIFEST.csv                                            #
    # ------------------------------------------------------------------ #
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    manifest_df = inv_df[["filename", "sha256", "bytes", "n_records", "source_acc"]].copy()
    manifest_df["created_utc"] = now_utc
    manifest_df["source_module"] = "00_raw_data"
    manifest_df["notes"] = ""

    # Reorder columns per INTERFACE.md
    manifest_df = manifest_df[
        ["filename", "sha256", "bytes", "n_records", "created_utc", "source_acc", "source_module", "notes"]
    ]

    out_path = MODULE_ROOT / "MANIFEST.csv"
    manifest_df.to_csv(out_path, index=False, lineterminator="\n")
    print(f"MANIFEST.csv written: {len(manifest_df)} rows → {out_path}")

    # ------------------------------------------------------------------ #
    # Step 5: Cross-check against phage_list.csv / bacteria_list.csv      #
    # 步骤 5：与 phage_list.csv / bacteria_list.csv 交叉核对                #
    # ------------------------------------------------------------------ #
    phage_list = pd.read_csv(MODULE_ROOT / "phage_list.csv")
    bacteria_list = pd.read_csv(MODULE_ROOT / "bacteria_list.csv")

    phage_dirs = {d.name for d in PHAGE_ROOT.iterdir() if d.is_dir()}
    bact_dirs = {d.name for d in BACTERIA_ROOT.iterdir() if d.is_dir()}

    issues = []

    missing_phage = set(phage_list["accession"]) - phage_dirs
    for acc in sorted(missing_phage):
        issues.append({"file": f"phage/{acc}/", "issue": "directory missing from disk", "recommended_fix": "Re-download via fetch_phages.py"})

    extra_phage = phage_dirs - set(phage_list["accession"])
    for acc in sorted(extra_phage):
        issues.append({"file": f"phage/{acc}", "issue": "file/dir not in phage_list.csv", "recommended_fix": "Investigate (may be pre-existing MANIFEST.csv)"})

    missing_bact = set(bacteria_list["accession"]) - bact_dirs
    known_invalid_bact = {"KY000037", "PY746849", "TODO"}
    for acc in sorted(missing_bact - known_invalid_bact):
        issues.append({"file": f"bacteria/{acc}/", "issue": "directory missing", "recommended_fix": "Check bacteria_list.csv notes"})

    # Check genome.fna completeness for phage dirs
    missing_genome_dirs = []
    for acc_dir in sorted(PHAGE_ROOT.iterdir()):
        if acc_dir.is_dir() and not (acc_dir / "genome.fna").exists():
            missing_genome_dirs.append(acc_dir.name)
    if missing_genome_dirs:
        for acc in missing_genome_dirs:
            issues.append({"file": f"phage/{acc}/genome.fna", "issue": "genome.fna missing", "recommended_fix": "Re-download"})

    if issues:
        issues_df = pd.DataFrame(issues)
        print(f"\nAudit issues ({len(issues_df)}):")
        print(issues_df.to_string(index=False))
    else:
        print("\nNo audit issues found.")

    return manifest_df, issues


if __name__ == "__main__":
    main()
