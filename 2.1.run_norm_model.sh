#!/bin/bash
# ============================================================================
# One array task per settings-combination; N worker processes per task run the
# CpGs (response) in parallel. Tailored to `midas` server:
#   - single node: 2x AMD EPYC 9654 = 192 physical cores / 384 threads, ~3 TB RAM
#   - partition `midas`, TIMELIMIT infinite
#   - QoS `normal`: MaxTRESPU cpu=80 -> 80 CPUs TOTAL across running jobs
#   - DefMemPerCPU=1024MB            -> always set --mem explicitly
#   - MaxArraySize=1001
#
# Sizing: 8 combos x 10 cpus = 80 = the per-user cap, so all 8 run at once.
#
#
# Environment (login node == compute node here, PyPI reachable):
#   cd ~/methylCHART && source .pcn_venv/bin/activate
#   pip install uv && uv pip install pcntoolkit   # optional: faster + lockfile
# ============================================================================
#SBATCH --job-name=methylCHART
#SBATCH --partition=midas
#SBATCH --array=0-3           # 4 combos (2 scales x 2 age forms)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8     # == number of worker processes
#SBATCH --threads-per-core=1  # 10 PHYSICAL cores: EPYC has 2 threads/core
#SBATCH --mem=32G             # DefMemPerCPU is only 1GB/cpu
#SBATCH --time=24:00:00       # partition allows `infinite`; kept as a guardrail
#SBATCH --output=logs/%x_%a.log
#SBATCH --error=logs/%x_%a.log

# set -euo pipefail
# mkdir -p logs

# One BLAS thread per worker: the parallelism is across CpGs, not in linear algebra. 
# Without this, N workers x M BLAS threads thrash the allocated cores.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONUNBUFFERED=1

cd "$HOME/methylCHART"
source .pcn_venv/bin/activate

# Keep all output on /home (local XFS)
# Note: savemodel=True produces many-small-files

# No srun: a single task keeps the full cpus-per-task allocation
python 2.1.run_norm_model.py "${SLURM_ARRAY_TASK_ID}"