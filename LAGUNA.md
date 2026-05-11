# LAGUNA.md — CARC Laguna HPC Cheat Sheet

For running heavy GPU jobs (ESM-2 3B embedding, AlphaFold 3 batch, Boltz-2 structure prediction). CPU-only smoke tests stay local.

**Cluster:** CARC Laguna (USC Center for Advanced Research Computing)
**OnDemand portal:** https://laguna-ood.carc.usc.edu (no VPN needed for CMC users)
**Account registration:** https://hpcaccount.usc.edu
**User:** CChen29@cmc.edu
**Project account (SLURM --account):** `jespinoza@kgi.edu_1541`
**Project directory:** `/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/`
**Conda env:** `/home1/CChen29@cmc.edu/.conda/envs/igem2026/` (permanent, persists across sessions)

---

## 1. Access via OnDemand (primary method — no VPN needed)

1. Go to **https://laguna-ood.carc.usc.edu**
2. Log in with CMC credentials (CChen29@cmc.edu)
3. Click **Code Server** in the left sidebar
4. Set parameters:
   - **Project account:** `jespinoza@kgi.edu_1541`
   - **Partition:** `compute` (for setup/job submission — does NOT need GPU)
   - **CPUs:** 2, **Memory:** 8 GB, **Hours:** 4
5. Click **Launch** → wait for session to start → click **Connect**
6. Inside VS Code: open terminal with **Ctrl+`**

**Every session: activate env and go to project directory**
```bash
source ~/.bashrc
conda activate igem2026
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
git pull
```

**Useful status checks:**
```bash
sinfo -p gpu                        # GPU partition availability (idle = ready)
squeue -j <jobid>                   # specific job status
sshare -A jespinoza@kgi.edu_1541    # allocation balance
sacct -u CChen29 --format=JobID,JobName,Elapsed,State   # job history
```

---

## 2. Submitting GPU jobs

> ⚠️ Always write the script to a file first, then submit.
> Do NOT paste heredoc (`<< 'EOF'`) with any leading spaces — the closing `EOF`
> must be at column 0 or it won't close.

**Standard pattern:**
```bash
sbatch scripts/boltz2_phiL7_tonB.slurm
```

Scripts are in `scripts/` and already configured. For custom jobs:
```bash
cat > /tmp/myjob.sh << 'EOF'
#!/bin/bash
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
...
EOF
sbatch /tmp/myjob.sh
```

**Monitor a running job:**
```bash
bash scripts/watch_job.sh <jobid>     # updates every 30s with log tail
tail -f logs/boltz2_<jobid>.out       # raw log
scancel <jobid>                       # cancel
```

---

## 3. Job templates

### Boltz-2 single pair — phiL7 tail spike × Xcc TonB (USE THIS FIRST)

Script already committed: `scripts/boltz2_phiL7_tonB.slurm`

```bash
git pull
sbatch scripts/boltz2_phiL7_tonB.slurm
```

The script automatically:
- Loads CUDA module (cuda-12.6.3)
- Installs torch>=2.7.0 from cu126 wheels (fixes `wrap_triton` issue)
- Falls back to trifast==0.1.10 if needed
- Runs Boltz-2 with 3 recycling steps, 200 sampling steps

Expected runtime: **15-30 min** on NVIDIA L40S. Output PDB:
```
05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/
```

### Boltz-2 batch screen (all RBP × receptor pairs)

```bash
cat > /tmp/boltz2_screen.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=boltz2_screen
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/logs/boltz2_screen_%j.out
#SBATCH --mail-user=CChen29@cmc.edu
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
mkdir -p logs
module load cuda 2>/dev/null
pip install "torch>=2.7.0" --index-url https://download.pytorch.org/whl/cu126 --force-reinstall -q

for fasta in 05_structure_prediction/inputs/boltz_input_*.fasta; do
    pair=$(basename "$fasta" .fasta | sed 's/boltz_input_//')
    boltz predict "$fasta" \
        --out_dir "05_structure_prediction/outputs/boltz2/$pair" \
        --accelerator gpu --recycling_steps 3 --sampling_steps 200 \
        --model boltz2 --use_msa_server --seed 42 --output_format pdb
done
EOF
sbatch /tmp/boltz2_screen.sh
```

### ESM-2 650M embedding (production quality)

```bash
cat > /tmp/esm2_650m.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=esm2_650m
#SBATCH --account=jespinoza@kgi.edu_1541
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/logs/esm2_%j.out
#SBATCH --mail-user=CChen29@cmc.edu
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate igem2026
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
mkdir -p logs

python 04_protein_embedding/processes/embed_esm2.py \
    --input 03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa \
    --model esm2_t33_650M_UR50D \
    --output 04_protein_embedding/outputs/embeddings_esm2_t33_650M_phiL7_rbps.npz \
    --batch-size 4 --pooling mean
EOF
sbatch /tmp/esm2_650m.sh
```

---

## 4. Pulling results back to local

**Recommended: OnDemand Files browser**
1. Go to https://laguna-ood.carc.usc.edu
2. Click **Files** in the top menu
3. Navigate to `/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/`
4. Select files → Download

**Alternative: rsync (once SSH is working)**
```bash
rsync -avz \
  "CChen29@laguna1.carc.usc.edu:/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/05_structure_prediction/outputs/" \
  05_structure_prediction/outputs/
```

---

## 5. Common gotchas

| Issue | Cause | Fix |
|-------|-------|-----|
| `CondaError: Run 'conda init'` | Shell not initialized | `source ~/.bashrc` first |
| `CUDA: False` in Code Server terminal | compute partition has no GPU | Expected — GPU only in sbatch jobs on `gpu` partition |
| Heredoc `EOF` not closing | Leading spaces before `EOF` | Use `cat > /tmp/job.sh << 'EOF'` with `EOF` at column 0 |
| `wrap_triton` ImportError | torch<2.7.0 — trifast 0.1.13 needs 2.7.0 | Use `scripts/boltz2_phiL7_tonB.slurm` which fixes this automatically |
| Job stuck in queue (`PD`) | GPU nodes busy | `sinfo -p gpu` — wait for idle nodes, or try `gpu_requeue` partition |
| `CUDA out of memory` | Batch too large | Reduce `--batch-size` to 2; or use ESM-2 150M |
| Permission denied on `/scratch/CChen29` | Wrong scratch path | Use `/project/jespinoza_1541/CChen29/` instead |

---

## 6. Storage

| Location | Purpose |
|----------|---------|
| `/home1/CChen29@cmc.edu/` | conda envs, `.bashrc` — persists forever |
| `/project/jespinoza_1541/CChen29/` | repo + outputs — use this for all work |

Keep raw genome data gitignored. Re-download on Laguna:
```bash
python 00_raw_data/processes/fetch_phages.py
python 00_raw_data/processes/fetch_bacteria.py
```

---

## 7. Budget tracking

After every GPU job, log to `08_cycle_data/outputs/hpc_log.csv`:
`date, job_name, jobid, gpu_hours, status, output_path, notes`

Check usage:
```bash
sacct -u CChen29 --format=JobID,JobName,Elapsed,State
```
