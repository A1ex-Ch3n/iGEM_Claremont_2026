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

## 1. Access (no VPN needed)

1. Go to https://laguna-ood.carc.usc.edu
2. Log in with CMC credentials (CChen29@cmc.edu)
3. Click **Clusters → Laguna Shell Access** for terminal
   OR click **Code Server** → Launch → open terminal with Ctrl+`

**SSH login node (once key is working):**
```bash
ssh CChen29@laguna1.carc.usc.edu
```

> SSH key fix: if portal didn't deploy key automatically, run in OnDemand terminal.
> Use `printf` (not echo) to avoid line-break issues when pasting:

```bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
printf 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIMEL3Fen+nBOJ8RHWK4ybEbAbGzXR2pSBEphodDU8CtU alex0071228@gmail.com\n' > ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
cat ~/.ssh/authorized_keys
```

**After opening terminal:**
```bash
source ~/.bashrc
conda activate igem2026
cd /project/jespinoza_1541/CChen29/iGEM_Claremont_2026
```

**Useful checks:**
```bash
sinfo -p gpu                       # see GPU nodes (partition name: gpu)
squeue -u CChen29                  # your jobs
sshare -A jespinoza@kgi.edu_1541   # allocation balance
```

---

## 2. Submitting jobs (write script to file first, then sbatch)

> ⚠️ Do NOT use heredoc (`<< 'EOF'`) with leading spaces — the closing EOF must be at
> column 0. Easiest fix: write the script to a file first, then submit.

**Pattern:**
```bash
cat > /tmp/myjob.sh << 'EOF'
#!/bin/bash
#SBATCH ...
...your commands...
EOF
sbatch /tmp/myjob.sh
```

**Monitor:**
```bash
squeue -u CChen29                                    # job status
tail -f logs/boltz2_<jobid>.out                      # live log
scancel <jobid>                                      # cancel a job
```

---

## 3. Job templates

### Template A — Boltz-2 single pair (phiL7 P25 × Xcc TonB)

```bash
cat > /tmp/boltz2_phiL7_tonB.sh << 'EOF'
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
EOF
sbatch /tmp/boltz2_phiL7_tonB.sh
```

Expected runtime: ~15-30 min on GPU. Output PDB at:
`05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`

### Template B — Boltz-2 batch screen (all RBP × receptor pairs)

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

for fasta in 05_structure_prediction/inputs/boltz_input_*.fasta; do
    pair=$(basename "$fasta" .fasta | sed 's/boltz_input_//')
    boltz predict "$fasta" \
        --out_dir "05_structure_prediction/outputs/boltz2/$pair" \
        --accelerator gpu \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --model boltz2 \
        --use_msa_server \
        --seed 42 \
        --output_format pdb
done
EOF
sbatch /tmp/boltz2_screen.sh
```

### Template C — ESM-2 650M embedding (production quality)

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
    --batch-size 4 \
    --pooling mean
EOF
sbatch /tmp/esm2_650m.sh
```

---

## 4. Pulling results back to local

After job completes, pull outputs from Laguna to your local machine.
Use the OnDemand Files browser, OR set up SSH key and rsync:

```bash
# Laguna → local: pull Boltz-2 structure results
rsync -avz \
  "CChen29@laguna1.carc.usc.edu:/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/05_structure_prediction/outputs/" \
  05_structure_prediction/outputs/

# Laguna → local: pull embedding results
rsync -avz \
  "CChen29@laguna1.carc.usc.edu:/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/04_protein_embedding/outputs/" \
  04_protein_embedding/outputs/
```

> SSH login hostname: check laguna-ood.carc.usc.edu → Help for the exact login node address.
> Alternatively use the OnDemand Files browser to download files directly.

---

## 5. Common gotchas

| Issue | Cause | Fix |
|-------|-------|-----|
| `CondaError: Run 'conda init'` | Shell not initialized | Run `source ~/.bashrc` first |
| `CUDA: False` in interactive session | compute partition has no GPU | Expected — GPU only available inside sbatch GPU jobs |
| Heredoc `EOF` not closing | Leading spaces before `EOF` | Write script to file with `cat > /tmp/job.sh << 'EOF'` then `sbatch /tmp/job.sh` |
| Job stuck in queue | GPU nodes busy | Check `sinfo -p gpu` for idle nodes; try `gpu_requeue` partition for preemptable jobs |
| `CUDA out of memory` | Batch too large | Drop `--batch-size` to 2; or use ESM-2 150M instead of 650M |
| Boltz-2 MSA timeout | MMSeqs2 server slow | Pre-computed MSAs already in `05_structure_prediction/outputs/boltz2/*/msa/` |

---

## 6. Project storage

| Location | Purpose | Size limit |
|----------|---------|-----------|
| `/home1/CChen29@cmc.edu/` | conda envs, personal configs | ~100 GB |
| `/project/jespinoza_1541/CChen29/` | repo + outputs (use this) | Shared with PI group |

Keep raw genome data OUT of the repo (gitignored). Re-download on Laguna with:
```bash
python 00_raw_data/processes/fetch_phages.py     # all 777 phages
python 00_raw_data/processes/fetch_bacteria.py   # all 34 bacteria
```

---

## 7. Budget tracking

After every GPU job, log to `08_cycle_data/outputs/hpc_log.csv`:
`date, job_name, jobid, gpu_hours, status, output_path, notes`

PI account `jespinoza@kgi.edu_1541` has limited SUs — don't waste them on test runs.
Use `sacct -u CChen29 --format=JobID,JobName,Elapsed,State` to review usage.
