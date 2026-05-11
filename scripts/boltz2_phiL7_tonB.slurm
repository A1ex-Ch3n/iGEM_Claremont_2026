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
nvidia-smi | grep -E "GPU Name|CUDA Version" | head -3

# Load CUDA runtime (confirmed: cuda-12.6.3 module available)
module load cuda 2>/dev/null
echo "[$(date)] CUDA_HOME=${CUDA_HOME:-NOT SET}"

# FIX (confirmed by local analysis):
# - torch 2.5.1+cu121: linux_x86_64 wheel, GCC-13 compatible, no CXXABI_1.3.15 issue
# - trifast==0.1.10: pure Python (py3-none-any), does NOT use wrap_triton
# This combination avoids both the wrap_triton and the CXXABI_1.3.15 errors.

echo "[$(date)] Installing torch 2.5.1+cu121..."
pip install "torch==2.5.1" \
    --index-url https://download.pytorch.org/whl/cu121 \
    --force-reinstall -q

echo "[$(date)] Downgrading trifast to 0.1.10 (no wrap_triton)..."
pip install "trifast==0.1.10" --force-reinstall -q

TORCH_VER=$(python -c "import torch; print(torch.__version__)")
CUDA_OK=$(python -c "import torch; print(torch.cuda.is_available())")
TRIFAST_OK=$(python -c "from trifast import triangle_attention; print('OK')" 2>&1)
echo "[$(date)] torch=$TORCH_VER  CUDA=$CUDA_OK  trifast=$TRIFAST_OK"

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
