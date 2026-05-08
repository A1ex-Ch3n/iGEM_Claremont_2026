#!/usr/bin/env bash
# Setup script: download PhageRBPdetect HMM database and models.
# 设置脚本：下载 PhageRBPdetect HMM 数据库和模型。
#
# Run from repo root: bash 03_rbp_identification/inputs/setup_inputs.sh
# 从仓库根目录运行：bash 03_rbp_identification/inputs/setup_inputs.sh
#
# Requirements: git, hmmpress (HMMER 3.4+)
# 要求：git, hmmpress（HMMER 3.4+）

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/phagerbpdetect_data"

echo "[1/4] Cloning PhageRBPdetection repository (shallow)..."
# 克隆 PhageRBPdetection 仓库（浅克隆）
git clone --depth 1 https://github.com/dimiboeckaerts/PhageRBPdetection.git /tmp/phagerbpdetect_src

echo "[2/4] Copying data files..."
# 复制数据文件
mkdir -p "$DATA_DIR"
cp /tmp/phagerbpdetect_src/data/RBPdetect_phageRBPs.hmm "$DATA_DIR/"
cp /tmp/phagerbpdetect_src/data/RBPdetect_xgb_model.json "$DATA_DIR/"
cp /tmp/phagerbpdetect_src/data/RBPdetect_xgb_hmm.json "$DATA_DIR/"
cp /tmp/phagerbpdetect_src/data/sequences.fasta "$DATA_DIR/"
rm -rf /tmp/phagerbpdetect_src

echo "[3/4] Pressing HMM database (required by hmmscan)..."
# Press HMM 数据库（hmmscan 所需）
HMMPRESS=$(which hmmpress 2>/dev/null || echo "/opt/homebrew/opt/hmmer/bin/hmmpress")
"$HMMPRESS" "$DATA_DIR/RBPdetect_phageRBPs.hmm"

echo "[4/4] Copying fallback phiL7 proteome..."
# 复制 phiL7 回退蛋白质组
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PHIL7_FAA="$REPO_ROOT/00_raw_data/phage/EU717894.1/proteins.faa"
if [ -f "$PHIL7_FAA" ]; then
    cp "$PHIL7_FAA" "$SCRIPT_DIR/EU717894.1_fallback_proteins.faa"
    echo "  Copied fallback FASTA from $PHIL7_FAA"
else
    echo "  WARNING: $PHIL7_FAA not found; run Module 00 first or provide manually."
fi

echo ""
echo "Done! Run: python3 03_rbp_identification/processes/run_pipeline.py"
# 完成！运行：python3 03_rbp_identification/processes/run_pipeline.py
