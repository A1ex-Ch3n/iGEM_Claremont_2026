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
conda activate boltz2
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
mkdir -p logs

echo "[$(date)] Job started on $(hostname)"
nvidia-smi | grep -E "GPU Name|CUDA Version" | head -3

module load cuda 2>/dev/null
echo "[$(date)] CUDA_HOME=${CUDA_HOME:-NOT SET}"

# Sanity check — skip trifast (it calls get_device_capability() at import,
# which fails outside GPU node; it will work fine when boltz calls it)
python -c "
import torch, optree
print(f'torch={torch.__version__}  CUDA={torch.cuda.is_available()}')
print(f'optree={optree.__version__}')
"

# Verify input sequence lengths before running
python -c "
seqs = {}; curr = None
for line in open('05_structure_prediction/inputs/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB.fasta'):
    if line.startswith('>'): curr = line.strip(); seqs[curr] = ''
    elif curr: seqs[curr] += line.strip()
for h, s in seqs.items():
    print(f'{h}: {len(s)} aa')
# Chain A must be 712 aa (rbp_01 tail spike), not 85 aa (P25)
"

# Clear the ENTIRE output directory (processed/ + msa/ + predictions/).
# Boltz-2 rebuilds processed/ from msa/ cache, so deleting only processed/
# still produces stale results. Full directory wipe forces fresh MSA download
# (~2 min) and fresh prediction.
OUT_DIR=05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB
rm -rf "${OUT_DIR}"
echo "[$(date)] Cleared all cached output. Running Boltz-2 from scratch..."

boltz predict \
  05_structure_prediction/inputs/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB.fasta \
  --out_dir "${OUT_DIR}" \
  --accelerator gpu \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --model boltz2 \
  --use_msa_server \
  --seed 42 \
  --output_format pdb

echo "[$(date)] Done. Exit code: $?"
