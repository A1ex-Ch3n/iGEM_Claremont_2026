#!/bin/bash
#SBATCH --job-name=boltz2_phiL7_tonB
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/logs/boltz2_%j.out
#SBATCH --mail-user=CChen29@cmc.edu
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
mkdir -p logs

echo "[$(date)] Job started on $(hostname)"
nvidia-smi | grep -E "NVIDIA-SMI|GPU Name|CUDA Version" | head -5

# Load CUDA runtime (CARC Laguna: cuda-12.6.3 module)
module load cuda 2>/dev/null
echo "[$(date)] CUDA_HOME=${CUDA_HOME:-NOT SET}"

# Root cause: wrap_triton requires torch>=2.7.0 (not 2.6.x).
# CARC has CUDA 12.6 so we use the cu126 PyTorch wheel index.
echo "[$(date)] Installing torch>=2.7.0 (cu126 wheels, compatible with CUDA 12.6)..."
pip install "torch>=2.7.0" \
    --index-url https://download.pytorch.org/whl/cu126 \
    --force-reinstall -q

TORCH_VER=$(python -c "import torch; print(torch.__version__)")
CUDA_OK=$(python -c "import torch; print(torch.cuda.is_available())")
echo "[$(date)] torch=$TORCH_VER  CUDA=$CUDA_OK"

# Verify wrap_triton exists (the actual root cause of previous failures)
python -c "from torch.library import wrap_triton; print('[$(date)] wrap_triton: OK')" || \
    echo "[$(date)] wrap_triton MISSING — trying Fix B: downgrade trifast to 0.1.10"
    pip install "trifast==0.1.10" --force-reinstall -q

echo "[$(date)] Running Boltz-2..."
boltz predict \
  05_structure_prediction/inputs/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB.fasta \
  --out_dir 05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB \
  --accelerator gpu \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --model boltz2 \
  --use_msa_server \
  --seed 42 \
  --output_format pdb

echo "[$(date)] Done. Exit code: $?"
