#!/bin/bash
# Run ONCE in Code Server terminal to create a clean boltz2 conda env.
# After this, all GPU jobs just activate boltz2 — no more env modification in jobs.
#
# Usage: bash scripts/setup_boltz2_env.sh
# Time: ~10 minutes

set -e  # exit on error

CONDA_BASE=/apps/conda/miniforge3/25.11.0-1
BOLTZ2_PREFIX=/home1/CChen29@cmc.edu/.conda/envs/boltz2
PIP=$BOLTZ2_PREFIX/bin/pip
PYTHON=$BOLTZ2_PREFIX/bin/python

echo "=== Creating fresh boltz2 conda env ==="
source $CONDA_BASE/bin/activate base
conda create -n boltz2 python=3.11 -y
echo "Done."

echo ""
echo "=== Installing packages in correct order ==="

# 1. optree FIRST — pin to 0.14.1 (manylinux_2_17, GCC-7 compat, no CXXABI_1.3.15)
#    Must be before torch so pip does NOT upgrade it
echo "[1/4] Installing optree==0.14.1..."
$PIP install "optree==0.14.1" -q

# 2. torch 2.5.1+cu121 — linux_x86_64 wheel, GCC-13 compat
#    optree is an 'extra' dep so pip will NOT upgrade it here
echo "[2/4] Installing torch 2.5.1+cu121..."
$PIP install "torch==2.5.1" \
    --index-url https://download.pytorch.org/whl/cu121 -q

# 3. boltz 2.0.3 — installs trifast>=0.1.11, we override in step 4
echo "[3/4] Installing boltz 2.0.3..."
$PIP install "boltz==2.0.3" -q

# 4. trifast 0.1.10 — pure Python (py3-none-any), does NOT use wrap_triton
echo "[4/4] Pinning trifast==0.1.10..."
$PIP install "trifast==0.1.10" --force-reinstall -q

echo ""
echo "=== Verification ==="
$PYTHON -c "
import torch; print('torch:', torch.__version__)
import optree; print('optree:', optree.__version__)
import trifast; print('trifast:', trifast.__version__)
import boltz; print('boltz:', boltz.__version__)
import pytorch_lightning; print('pytorch_lightning:', pytorch_lightning.__version__)
print()
print('All imports OK.')
print('NOTE: CUDA=False here (login node has no GPU).')
print('CUDA will be True on GPU node when job runs.')
"

echo ""
echo "=== Setup complete ==="
echo "Now submit jobs with: sbatch scripts/boltz2_phiL7_tonB.slurm"
