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
nvidia-smi | head -12

# Load CUDA runtime so PyTorch can find libcuda.so
echo "[$(date)] Loading CUDA module..."
module load cuda 2>/dev/null || \
module load cuda/12.4 2>/dev/null || \
module load cuda/12.1 2>/dev/null || \
echo "WARNING: no cuda module found — trying anyway"

echo "[$(date)] CUDA_HOME=$CUDA_HOME  LD_LIBRARY_PATH (first 3): $(echo $LD_LIBRARY_PATH | tr ':' '\n' | head -3 | tr '\n' ':')"

# trifast 0.1.13 requires torch>=2.6.0 (wrap_triton added in 2.6)
echo "[$(date)] Installing GPU PyTorch >=2.6.0..."
pip install "torch>=2.6.0" --index-url https://download.pytorch.org/whl/cu121 -q --force-reinstall

echo "[$(date)] CUDA available: $(python -c 'import torch; print(torch.cuda.is_available())')"
echo "[$(date)] GPU count: $(python -c 'import torch; print(torch.cuda.device_count())')"

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
