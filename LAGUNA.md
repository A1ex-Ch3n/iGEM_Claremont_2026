# LAGUNA.md — HPC Operation Cheat Sheet

For running heavy GPU jobs (ESM-2 3B embedding, AlphaFold 3 batch, Boltz-2 large screens) on **CARC Laguna** (USC Center for Advanced Research Computing regional cluster). CPU-only smoke tests stay local.

**Cluster:** CARC Laguna  
**OnDemand portal:** https://laguna-ood.carc.usc.edu (browser-based, no VPN needed)  
**Account registration:** https://hpcaccount.usc.edu  
**User:** CChen29@cmc.edu  
**Project account (SLURM --account):** `jespinoza@kgi.edu_1541`  
**SSH key:** deployed 2026-05-10

---

## 1. Access

**Recommended (no VPN needed):** Use the OnDemand web portal

1. Go to https://laguna-ood.carc.usc.edu
2. Log in with CMC credentials
3. Click **Clusters → Laguna Shell Access** for a browser terminal
4. Or click **JupyterLab** in the left sidebar for a notebook interface

**SSH (requires being on campus or USC VPN):**
```bash
# Exact hostname TBD — check laguna-ood.carc.usc.edu → Help for SSH details
ssh CChen29@<laguna-login-node>
```

**First-time checks in terminal:**
```bash
sinfo -p gpu                          # see GPU partition availability
sshare -A jespinoza@kgi.edu_1541      # see PI's allocation balance
squeue -u CChen29                     # see your jobs
```

**Useful links:**
- User guide: https://www.carc.usc.edu/user-guides
- Storage: `$SCRATCH` = `/scratch/CChen29` (large, not backed up)

---

## 2. Environment setup on Laguna (one-time)

```bash
# Pick a working directory under your scratch / project space
mkdir -p $SCRATCH/igem_2026 && cd $SCRATCH/igem_2026

# Clone the repo
git clone <your-fork-or-team-remote> .
git checkout active-learning-pipeline

# Load module(s) — names are cluster-specific, verify with `module avail`
module load anaconda3
module load cuda/12.1     # or whatever current version

# Create env from the project file
conda env create -f shared/env/environment.yml
conda activate igem2026

# GPU-only extras (CUDA build of PyTorch + ESM)
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install fair-esm
```

Test GPU is visible:

```bash
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.device_count())"
```

---

## 3. SLURM job templates

### Template A — single-GPU ESM-2 3B embedding

`scripts/slurm/embed_esm2_3b.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=esm2_3b
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu           # verify with: sinfo -p gpu           # e.g., gpu, gpu-a100
#SBATCH --gres=gpu:a100:1                     # adjust to available GPU
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH [email protected]
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
module load cuda/12.1

cd $SCRATCH/igem_2026

python 04_protein_embedding/processes/embed_esm2.py \
    --input 03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa \
    --model esm2_t36_3B_UR50D \
    --output 04_protein_embedding/outputs/embeddings_esm2_t36_3B_phiL7_rbps.npz \
    --batch-size 4 \
    --pooling mean
```

Submit: `sbatch scripts/slurm/embed_esm2_3b.slurm`

Monitor: `squeue -u $USER` · `tail -f logs/esm2_3b_<jobid>.out`

### Template B — multi-GPU AlphaFold 3 batch

`scripts/slurm/af3_batch.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=af3_batch
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu           # verify with: sinfo -p gpu
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH [email protected]
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
module load cuda/12.1

cd $SCRATCH/igem_2026

# AF3 setup is non-trivial — see the per-module SETUP_BLOCKERS.md
# for the current install state. Fallback: ColabFold via local API.
python 05_structure_prediction/processes/run_af3.py \
    --input 03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa \
    --output 05_structure_prediction/outputs/af3/ \
    --max-templates 4 \
    --num-recycles 3
```

### Template C — Boltz-2 affinity screen

`scripts/slurm/boltz2_screen.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=boltz2_screen
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu           # verify with: sinfo -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH [email protected]
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
module load cuda/12.1

cd $SCRATCH/igem_2026

python 05_structure_prediction/processes/run_boltz2.py \
    --rbps 03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa \
    --receptors data/xcc_receptors.faa \
    --output 05_structure_prediction/outputs/boltz2/ \
    --predict-affinity
```

---

## 4. File transfer (local ⇄ Laguna)

```bash
# Local → Laguna: push the latest code
rsync -avz --exclude='00_raw_data/phage' --exclude='00_raw_data/bacteria' \
    "$(pwd)/" "CChen29@<laguna-login-node>:$SCRATCH/igem_2026/"

# Laguna → local: pull computed outputs
rsync -avz "CChen29@<laguna-login-node>:$SCRATCH/igem_2026/04_protein_embedding/outputs/" \
    04_protein_embedding/outputs/
```

Don't `rsync` the raw genome data; it's reproducible from NCBI.

---

## 5. Common gotchas

| Issue | Cause | Fix |
|---|---|---|
| `CUDA out of memory` | Batch too large for ESM-2 3B on A100 | Drop `--batch-size` to 2; or use ESM-2 650M |
| Job stuck in queue | Wrong partition name | `sinfo -s` to list partitions |
| `module: command not found` | Login shell doesn't source modules | `source /etc/profile.d/modules.sh` first |
| AF3 install fails | Model weight access not granted | Use ColabFold (`pip install colabfold`) as drop-in fallback for tonight |
| Boltz-2 install fails | Missing nvcc | Load CUDA module before pip install |

---

## 6. Notification setup

Every SLURM job above includes `#SBATCH [email protected]` and `--mail-type=END,FAIL`. You'll get email when jobs finish or crash.

---

## 7. Project budget tracking

After every big run, log to `08_cycle_data/outputs/hpc_log.csv` with columns:
`date, job_name, jobid, gpu_hours, status, output_path, notes`

So we can stay within the academic allocation (~200 GPU-hours estimated for project duration).
